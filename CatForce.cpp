// CatForce - Catalyst search utility based on LifeAPI using brute force.
// Written by Michael Simkin 2015
#include "LifeAPI.h"
#include <algorithm>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>

const int MAX_CATALYSTS = 5;

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
}

class SearchParams {
public:
  int maxGen;
  //int numCatalysts; // TODO: remove?
  int numTransparent;
  int stableInterval;
  std::string activeRegion;
  int xActiveRegion;
  int yActiveRegion;
  int startGen;
  int lastGen;
  std::string outputFile;
  std::string fullReportFile;
  int searchArea[4]{};
  int maxW;
  int maxH;
  std::vector<std::string> targetFilter;
  std::vector<int> filterdx;
  std::vector<int> filterdy;
  std::vector<int> filterGen;

  //modification
  std::vector<StaticSymmetry> filterGroups;
  // number of generations: if catalysts are destroyed, filters must be met within this many generations afterward or earlier.
  int stopAfterCatsDestroyed; 
  bool reportMatches;
  bool quietMode;

  std::vector<std::pair<int, int>> filterGenRange;
  //std::tuple<std::string, int, int> alsoRequired;

  int maxCatSize;


  std::vector<int> offsetBboxes;
  std::pair<int,int> objReoccurs;
  std::pair<int,int> numCatalystsPrePost;

  int startSymInteraction;
  int lastSymInteraction;

  SearchParams() {
    maxGen = 250;
    // numCatalysts = 2;
    stableInterval = 15;
    activeRegion = "";
    searchArea[0] = -10;
    searchArea[1] = 0;
    searchArea[2] = 20;
    searchArea[3] = 20;
    xActiveRegion = 0;
    yActiveRegion = 0;
    startGen = 1;
    lastGen = -1;
    outputFile = "results.rle";
    maxW = -1;
    maxH = -1;
    //symmetryEnum = StaticSymmetry::C1;
    maxCatSize = -1;
    fullReportFile = "";

    //alsoRequired = std::make_tuple("", 0,0);
    
    filterGroups = {};
    stopAfterCatsDestroyed = -1;

    quietMode = false;
    reportMatches = false;

    startSymInteraction = 5;
    lastSymInteraction = 25;
    offsetBboxes = {};
    numCatalystsPrePost = std::make_pair(2,1);
    objReoccurs = std::make_pair(0,0);
  }
};

class CatalystInput {
public:
  std::string rle;
  int maxDisappear;
  int centerX;
  int centerY;
  char symmType;
  std::vector<std::string> forbiddenRLE;
  std::vector<std::pair<int, int>> forbiddenXY;
  std::string requiredRLE;
  std::pair<int, int> requiredXY;
  std::string locusRLE;
  std::pair<int, int> locusXY;
  bool transparent;

  explicit CatalystInput(std::string &line) {
    std::vector<std::string> elems;
    split(line, ' ', elems);

    if (elems.size() < 6) {
      std::cout << "The line " << line << "is invalid" << std::endl;
      std::cout << "Format: cat <rle> <absense interval> <centerX> <centerY> "
                   "<symm Type | + / x *>"
                << std::endl;
      exit(1);
    }

    rle = elems[1];
    maxDisappear = atoi(elems[2].c_str());
    centerX = atoi(elems[3].c_str());
    centerY = atoi(elems[4].c_str());
    symmType = elems[5].at(0);

    transparent = false;

    int argi = 6;

    while (argi < elems.size()) {
      if (elems[argi] == "forbidden") {
        forbiddenRLE.push_back(elems[argi + 1]);
        forbiddenXY.emplace_back(
            atoi(elems[argi + 2].c_str()), atoi(elems[argi + 3].c_str()));

        argi += 4;
      } else if (elems[argi] == "required") {
        requiredRLE = elems[argi + 1];
        requiredXY = std::make_pair(atoi(elems[argi + 2].c_str()), atoi(elems[argi + 3].c_str()));
        argi += 4;
      } else if (elems[argi] == "locus") {
        locusRLE = elems[argi + 1];
        locusXY = std::make_pair(atoi(elems[argi + 2].c_str()), atoi(elems[argi + 3].c_str()));
        argi += 4;
      } else if (elems[argi] == "transparent") {
        transparent = true;
        argi += 1;
      } else {
        std::cout << "Unknown catalyst attribute: " << elems[argi] << std::endl;
        exit(1);
      }
    }
  }

  void Print() const {
    std::cout << rle << " " << maxDisappear << " " << centerX << " " << centerY
              << " " << symmType << std::endl;
  }
};

void CharToTransVec(char ch, std::vector<AffineTransform> &trans) {

  if (ch == '.'){
    trans = SymmetryGroupFromEnum(StaticSymmetry::C1);
    return;
  } if (ch == '|') {
    trans =  SymmetryGroupFromEnum(StaticSymmetry::D2AcrossY);
    return;
  }

  if (ch == '-') {
    trans =  SymmetryGroupFromEnum(StaticSymmetry::D2AcrossX);
    return;
  }

  if (ch == '+') {
    trans =  SymmetryGroupFromEnum(StaticSymmetry::C4);
    return;
  }

  if (ch == '/') {
    trans =  SymmetryGroupFromEnum(StaticSymmetry::D2negdiagodd);
    return;
  }
  // For 180 degree symetrical
  if (ch == 'x') {
    trans.push_back(AffineTransform());
    trans.push_back(AffineTransform(LinearTransform::Rotate90));
    trans.push_back(AffineTransform(LinearTransform::FlipAcrossX));
    trans.push_back(AffineTransform(LinearTransform::FlipAcrossYEqX));
    return;
  }

  if (ch == '*') {
    trans=SymmetryGroupFromEnum(StaticSymmetry::D8);
    return;
  }
}


void ReadParams(const std::string& fname, std::vector<CatalystInput> &catalysts,
                SearchParams &params) {
  std::ifstream infile;
  infile.open(fname.c_str(), std::ifstream::in);
  if (!infile.good()) {
    std::cout << "Could not open file!" << std::endl;
    exit(1);
  }

  std::string Cat = "cat";
  std::string maxGen = "max-gen";
  std::string startGen = "start-gen";
  std::string lastGen = "last-gen";

  std::string numTransp = "num-transparent";
  std::string stable = "stable-interval";
  std::string area = "search-area";
  std::string outputFile = "output";
  std::string filter = "filter";
  std::string maxWH = "fit-in-width-height";
  std::string maxCatSize = "max-category-size";
  std::string fullReport = "full-report";
  std::string alsoRequired = "also-required";

  std::string stopAfterCatsDestroyed = "stop-after-cats-destroyed";
  std::string reportMatches = "report-matches";
  std::string quiet = "quiet-mode";

  std::string line;
  std::string activeReg = "active-region";
  std::string objReoccurs = "object-reoccurs";
  std::string numCatPreSym = "num-cats-pre-sym";
  std::string numCatPostSym = "num-cats-post-sym";
  std::string startSymmetricInteraction = "start-sym-interaction";
  std::string lastSymmetricInteraction = "last-sym-interaction";
  std::string offsets = "offsets";

  while (std::getline(infile, line)) {
    try {
      std::vector<std::string> elems;
      split(line, ' ', elems);

      if (elems.size() < 2)
        continue;

      if (elems[0] == Cat)
        catalysts.emplace_back(line);

      if (elems[0] == maxGen)
        params.maxGen = atoi(elems[1].c_str());

      if (elems[0] == numTransp)
        params.numTransparent = atoi(elems[1].c_str());

      if (elems[0] == stable)
        params.stableInterval = atoi(elems[1].c_str());

      if (elems[0] == activeReg) {
        params.activeRegion = elems[1];

        if (elems.size() > 3) {
          params.xActiveRegion = atoi(elems[2].c_str());
          params.yActiveRegion = atoi(elems[3].c_str());
        }
      }

      if (elems[0] == area) {
        params.searchArea[0] = atoi(elems[1].c_str());
        params.searchArea[1] = atoi(elems[2].c_str());
        params.searchArea[2] = atoi(elems[3].c_str());
        params.searchArea[3] = atoi(elems[4].c_str());
      }

      if (elems[0] == offsets){
        assert(elems.size() >= 5);
        params.offsetBboxes.push_back(atoi(elems[1].c_str()));
        params.offsetBboxes.push_back(atoi(elems[2].c_str()));
        params.offsetBboxes.push_back(atoi(elems[3].c_str()));
        params.offsetBboxes.push_back(atoi(elems[4].c_str()));
      }

      if (elems[0] == numCatPreSym ){
        params.numCatalystsPrePost.first = atoi(elems[1].c_str());
      }

      if (elems[0] == numCatPostSym ){
        params.numCatalystsPrePost.second = atoi(elems[1].c_str());
      }

      if (elems[0] == startGen)
        params.startGen = atoi(elems[1].c_str());

      if (elems[0] == lastGen)
        params.lastGen = atoi(elems[1].c_str());

      if (elems[0] == outputFile) {
        params.outputFile = elems[1];

        for (int i = 2; i < elems.size(); i++) {
          params.outputFile.append(" ");
          params.outputFile.append(elems[i]);
        }
      }

      if (elems[0] == fullReport) {
        params.fullReportFile = elems[1];

        for (int i = 2; i < elems.size(); i++) {
          params.fullReportFile.append(" ");
          params.fullReportFile.append(elems[i]);
        }
      }

      if (elems[0] == startSymmetricInteraction)
        params.startSymInteraction = atoi(elems[1].c_str());

      if (elems[0] == lastSymmetricInteraction)
        params.lastSymInteraction = atoi(elems[1].c_str());

      if (elems[0] == objReoccurs){
        std::vector<std::string> rangeElems;
        split(elems[1], '-', rangeElems);
        if (rangeElems.size() == 1 && elems.size() < 3){
          std::cout << "input for object-reoccurs not understood" << std::endl;
          exit(0);
        } else if (elems.size() >= 3) {
          params.objReoccurs = std::make_pair(atoi(elems[1].c_str()), atoi(elems[2].c_str()));
        } else {
          params.objReoccurs = std::make_pair(atoi(rangeElems[0].c_str()), atoi(rangeElems[1].c_str()));
        }
      }

      if (elems[0] == filter) {
        std::vector<std::string> rangeElems;
        split(elems[1], '-', rangeElems);

        if (rangeElems.size() == 1) {
          params.filterGen.push_back(atoi(elems[1].c_str()));
          params.filterGenRange.emplace_back(-1, -1);
        } else {
          int minGen = atoi(rangeElems[0].c_str());
          int maxGen = atoi(rangeElems[1].c_str());

          params.filterGen.push_back(-1);
          params.filterGenRange.emplace_back(minGen, maxGen);
        }

        params.targetFilter.push_back(elems[2]);
        params.filterdx.push_back(atoi(elems[3].c_str()));
        params.filterdy.push_back(atoi(elems[4].c_str()));

        // modification
        if ( elems.size() >= 6 && (elems[5].at(0) == 'D' || elems[5].at(0) == 'C') ) {
          if ( elems[5] == "C1" || SymmetryEnumFromString(elems[5]) != StaticSymmetry::C1){
            params.filterGroups.push_back(SymmetryEnumFromString(elems[5]));
          } else {
            std::cout << "filter symmetry group " << elems[5] << " not understood." << std::endl;
            exit(0);
          }
        } else {
          params.filterGroups.push_back(StaticSymmetry::C1);
        }

      }

      if (elems[0] == maxWH) {
        params.maxW = atoi(elems[1].c_str());
        params.maxH = atoi(elems[2].c_str());
      }

      if (elems[0] == maxCatSize) {
        params.maxCatSize = atoi(elems[1].c_str());
      }

      if (elems[0] == stopAfterCatsDestroyed){
        params.stopAfterCatsDestroyed = atoi(elems[1].c_str());
      }

      if (elems[0] == reportMatches){
        params.reportMatches = true;
      }




      if (elems[0] == quiet){
        params.quietMode = true;
      }

    } catch (const std::exception &ex) {
    }
  }
  if(params.lastGen == -1)
    params.lastGen = params.maxGen - 1;

  if (params.activeRegion.length() == 0) {
    std::cout << "Did not read any pattern!" << std::endl;
    exit(1);
  }
  if (catalysts.empty()) {
    std::cout << "Did not read any catalysts!" << std::endl;
    exit(1);
  }
}

struct Configuration {
  int count;
  std::pair<int,int> prePostCount;
  int transparentCount;
  std::array<int, MAX_CATALYSTS> curx;
  std::array<int, MAX_CATALYSTS> cury;
  std::array<int, MAX_CATALYSTS> curs;
  // int minIter;
  LifeState state;
  LifeState catalystsState;
  // debugging purposes
  LifeState catsPostSymmetry;
  LifeState catsPreSymmetry;
  // end debugging purposes
  bool postSymmetry;
  std::pair<int,int> loneOffset;
  // std::vector<LifeTarget> shiftedTargets;
};

// Fix a, what positions of b causes a collision?
LifeState CollisionMask(const LifeState &a, const LifeState &b) {
  int popsum = a.GetPop() + b.GetPop();

  LifeState mask;
  for (int x = 0; x < N; x++) {
    for (int y = 0; y < 64; y++) {
      LifeState state = a;
      state.Join(b, x, y);

      // No overlaps allowed
      if (!(state.GetPop() == popsum)) {
        mask.Set(x, y);
        continue;
      }

      state.Step();

      if (!state.Contains(a) || !state.Contains(b, x, y)) {
        mask.Set(x, y);
      }

    }
  }

  return mask;
}

std::string GetRLE(const LifeState &s);

LifeState LoadCollisionMask(const LifeState &a, const LifeState &b) {
  std::stringstream ss;
  ss << "masks/mask-" << a.GetHash() << "-" << b.GetHash();
  std::string fname = ss.str();

  std::ifstream infile;
  infile.open(fname.c_str(), std::ifstream::in);
  if (!infile.good()) {
    LifeState mask = CollisionMask(a, b);
    std::ofstream outfile;
    outfile.open(fname.c_str(), std::ofstream::out);
    outfile << GetRLE(mask);
    return mask;
  } else {
    std::stringstream buffer;
    buffer << infile.rdbuf();
    std::string rle = buffer.str();
    LifeState result = LifeState::Parse(rle.c_str());
    result.Move(-32, -32);
    return result;
  }
}

void GenerateStates(const std::vector<CatalystInput> &catalysts,
                    std::vector<LifeState> &states,
                    std::vector<LifeState> &required,
                    std::vector<LifeState> &locus,
                    std::vector<std::vector<LifeTarget>> &forbidden,
                    std::vector<int> &maxMissing,
                    std::vector<bool> &transparent) {

  for (const auto & catalyst : catalysts) {
    std::vector<AffineTransform> trans;
    CharToTransVec(catalyst.symmType, trans);

    const char *rle = catalyst.rle.c_str();
    int dx = catalyst.centerX;
    int dy = catalyst.centerY;
    int maxDisappear = catalyst.maxDisappear;

    for (auto & tran : trans) {
      LifeState pat = LifeState::Parse(rle, dx, dy, tran);
      states.push_back(pat);
      maxMissing.push_back(maxDisappear);

      std::vector<LifeTarget> forbidTarg;

      for (int k = 0; k < catalyst.forbiddenRLE.size(); k++) {
        forbidTarg.push_back(LifeTarget::Parse(catalyst.forbiddenRLE[k].c_str(),
                                               catalyst.forbiddenXY[k].first,
                                               catalyst.forbiddenXY[k].second, tran));
      }

      forbidden.push_back(forbidTarg);

      LifeState catrequired;
      if (catalyst.requiredRLE != "") {
        catrequired = LifeState::Parse(catalyst.requiredRLE.c_str(),
                                       catalyst.requiredXY.first,
                                       catalyst.requiredXY.second, tran);
      }
      required.push_back(catrequired);

      LifeState catlocus;
      if (catalyst.locusRLE != "") {
        catlocus = LifeState::Parse(catalyst.locusRLE.c_str(),
                                       catalyst.locusXY.first,
                                       catalyst.locusXY.second, tran);
      } else {
        catlocus = pat;
      }
      locus.push_back(catlocus);

      transparent.push_back(catalyst.transparent);
    }
  }
}

std::string GetRLE(const std::vector<std::vector<int>> &life2d) {
  if (life2d.empty())
    return "";

  if (life2d[0].empty())
    return "";

  int h = life2d[0].size();

  std::stringstream result;

  int eol_count = 0;

  for (int j = 0; j < h; j++) {
    int last_val = -1;
    int run_count = 0;

    for (const auto & i : life2d) {
      int val = i[j];

      // Flush linefeeds if we find a live cell
      if (val == 1 && eol_count > 0) {
        if (eol_count > 1)
          result << eol_count;

        result << "$";

        eol_count = 0;
      }

      // Flush current run if val changes
      if (val == 1 - last_val) {
        if (run_count > 1)
          result << run_count;

        if (last_val == 1)
          result << "o";
        else
          result << "b";

        run_count = 0;
      }

      run_count++;
      last_val = val;
    }

    // Flush run of live cells at end of line
    if (last_val == 1) {
      if (run_count > 1)
        result << run_count;

      result << "o";

      run_count = 0;
    }

    eol_count++;
  }

  // Flush trailing linefeeds
  if (eol_count > 0) {
    if (eol_count > 1)
      result << eol_count;

    result << "$";

    eol_count = 0;
  }

  return result.str();
}

std::string GetRLE(const LifeState &s) {
  std::vector<std::vector<int>> vec(N, std::vector<int>(N));

  for (int j = 0; j < N; j++)
    for (int i = 0; i < N; i++)
      vec[i][j] = s.GetCell(i - 32, j - 32);

  return GetRLE(vec);
}

class SearchResult {
public:
  // Saved for the report
  LifeState init;

  // iters state in form of integers
  // std::vector<int> params;
  int maxGenSurvive;
  int firstGenSurvive;

  SearchResult(LifeState &initState, const Configuration &conf,
               int firstGenSurviveIn, int genSurvive) {
    init.Copy(initState);
    maxGenSurvive = genSurvive;
    firstGenSurvive = firstGenSurviveIn;
  }

  void Print() {
    std::cout << "start:" << firstGenSurvive;
    std::cout << ", finish:" << maxGenSurvive << ", params: ";
    std::cout << std::endl;
  }
};

// TODO: change this to group by generation of results [and offset too?]
// and sort by amount of junk left over in that generation.
class Category {
private:
  int catDelta;
  int maxgen;
  uint64_t hash;

public:
  LifeState categoryKey;
  std::vector<SearchResult> results;

  Category(LifeState &catalystRemoved, SearchResult &firstResult,
           int catDeltaIn, int maxGen) {
    categoryKey.Copy(catalystRemoved);
    results.push_back(firstResult);
    catDelta = catDeltaIn;
    maxgen = maxGen;

    LifeState temp = categoryKey;
    temp.Step(maxgen - temp.gen);
    hash = temp.GetHash();
  }

  void Add(SearchResult &result) { results.push_back(result); }

  bool BelongsTo(LifeState &test, const uint64_t &testHash) {
    if (testHash != hash)
      return false;

    LifeState tempCat = categoryKey;
    LifeState tempTest = test;

    if (tempTest.gen > tempCat.gen)
      tempCat.Step(tempTest.gen - tempCat.gen);
    else if (tempTest.gen < tempCat.gen)
      tempTest.Step(tempCat.gen - tempTest.gen);

    for (int i = 0; i < catDelta; i++) {
      if (tempTest == tempCat)
        return true;

      tempCat.Step();
      tempTest.Step();
    }

    return false;
  }

  static bool CompareSearchResult(SearchResult &a, SearchResult &b) {
    return (a.maxGenSurvive - a.firstGenSurvive) >
           (b.maxGenSurvive - b.firstGenSurvive);
  }

  void Sort() {
    std::sort(results.begin(), results.end(), Category::CompareSearchResult);
  }

  void Print() {
    for (auto & result : results)
      result.Print();
  }

  std::string RLE(int maxCatSize) {
    // 36 is extra margin to get 100
    const int Dist = 36 + 64;

    int width = Dist * results.size();
    int height = Dist;

    std::vector<std::vector<int>> vec;

    for (int i = 0; i < width; i++) {
      std::vector<int> temp;

      for (int j = 0; j < height; j++)
        temp.push_back(0);

      vec.push_back(temp);
    }

    int howmany = results.size();

    if (maxCatSize != -1)
      howmany = std::min(howmany, maxCatSize);

    for (int l = 0; l < howmany; l++)
      for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++)
          vec[Dist * l + i][j] = results[l].init.GetCell(i - 32, j - 32);

    return GetRLE(vec);
  }
};

class CategoryContainer {
public:
  std::vector<Category *> categories;
  int catDelta;
  int maxgen;

  explicit CategoryContainer(int maxGen) {
    catDelta = 14;
    maxgen = maxGen + catDelta;
  }

  CategoryContainer(int cats, int maxGen) {
    catDelta = cats;
    maxgen = maxGen + catDelta;
  }

  void Add(LifeState &init, LifeState &afterCatalyst, LifeState &catalysts,
           const Configuration &conf, int firstGenSurvive,
           int genSurvive) {
    LifeState result;

    result.Copy(afterCatalyst);
    result.Copy(catalysts, XOR);

    result.Step(maxgen - result.gen);
    uint64_t hash = result.GetHash(); // we group by active cells at maxGen

    result.Clear();
    result.Copy(afterCatalyst);
    result.Copy(catalysts, XOR);

    // see if it fits in an existing category.
    for (auto & category: categories) {
      if (category->BelongsTo(result, hash)) {
          SearchResult r(init, conf, firstGenSurvive, genSurvive);
          category->Add(r);
          return;
      }
    }

    // create new category.
    LifeState categoryKey;
    categoryKey.Copy(afterCatalyst);
    categoryKey.Copy(catalysts, XOR);

    SearchResult r(init, conf, firstGenSurvive, genSurvive);
    categories.push_back(new Category(categoryKey, r, catDelta, maxgen));
  }

  void Sort() {
    for (auto & category: categories)
      category->Sort();
  }

  void Print() {
    for (auto & category: categories)
      category->Print();
  }

  std::string CategoriesRLE(int maxCatSize) {
    std::stringstream ss;
    for (auto & category: categories) {
      ss << category->RLE(maxCatSize);
    }

    return ss.str();
  }
};


class CatalystSearcher {
public:
  clock_t begin{};

  SearchParams params;

  LifeState activeRegion;
  LifeState offsets;
  LifeState searchArea;
  int maxNumCatalysts{};
  
  std::vector<LifeState> catalysts;
  std::vector<int> maxMissing;
  std::vector<LifeTarget> targets;
  std::vector<std::vector<LifeTarget>> forbiddenTargets;
  std::vector<LifeState> requiredParts;
  std::vector<std::vector<LifeState>> catalystCollisionMasks;
  std::vector<LifeState> catalystReactionMasks;
  std::vector<LifeState> catalystLocus;
  std::vector<LifeState> catalystLocusReactionMasks;
  std::vector<LifeState> catalystAvoidMasks;
  std::vector<bool>      transparent;

  std::vector<std::array<int,3>> rotatedCatalystMatches;
  std::vector<LinearTransform> catTransformations;
  
  clock_t current{};
  time_t lastReport{};
  long long idx{};
  int found{};
  int fullfound{};
  long long total{};
  unsigned short int counter{};
  
  CategoryContainer *categoryContainer{};
  CategoryContainer *fullCategoryContainer{};

  // flags and memeber for the search
  
  //std::vector<LifeTarget> targetFilter;
  std::vector<std::vector<LifeTarget>> targetFilterLists;
  std::vector<int> matches;
  bool hasFilter{};
  bool reportAll{};
  bool hasFilterDontReportAll{};
  int filterMaxGen{};



  void GenerateOffsets(){
    for(int i = 0; i < params.offsetBboxes.size(); i+= 4){
      offsets.Join(LifeState::SolidRect(params.offsetBboxes[i], params.offsetBboxes[i+1],
                              params.offsetBboxes[i+2], params.offsetBboxes[i+3]));
    }
    // eliminate those such that active region overlaps with its image
    LifeTarget activeRegionTarget(activeRegion);
    std::array<int,4> xyBounds = activeRegion.XYBounds();
    for(int x = 2*xyBounds[0]-2; x <= 2*xyBounds[2]+2; ++x){
      for(int y = 2*xyBounds[1]-2; y <= 2*xyBounds[3]+2; ++y){
        if (offsets.GetCell(x,y) == 1){
          LifeState combined = activeRegion;
          combined.Transform(Rotate180OddBoth);
          combined.Move(x,y);
          combined.Join(activeRegion);
          if(!combined.Contains(activeRegionTarget)){
            offsets.Erase(x,y);
          }
        }
      }
    }
  }

  void ComputeCatalystTransformationData(std::vector<CatalystInput> & catInput){
    for(auto catData : catInput){
      std::vector<AffineTransform> whichTransforms;
      CharToTransVec(catData.symmType, whichTransforms);
      for (AffineTransform T : whichTransforms){
        catTransformations.push_back(T.matrix);
      }
    }
    int unorientedCatIndex = -1;
    int defaultCatIndex = -1;
    for(int catIndex = 0; catIndex < catalysts.size(); ++catIndex){
      if(catTransformations[catIndex] == Rotate0){
        defaultCatIndex = catIndex;
        unorientedCatIndex += 1;
      }
      LifeState transformedCatalyst = catalysts[catIndex];
      transformedCatalyst.Transform(Rotate180OddBoth);
      for( int i = defaultCatIndex; i < std::min(size_t(defaultCatIndex+8), catalysts.size()); ++i){
        LifeState candidateMatch = catalysts[i];
        candidateMatch.Move(32,32);
        transformedCatalyst.Move(32,32);
        std::pair<int,int> transformedFirstOn = transformedCatalyst.FirstOn();
        std::pair<int,int> candidateFirstOn = candidateMatch.FirstOn();
        int xDiff = transformedFirstOn.first-candidateFirstOn.first;
        int yDiff = transformedFirstOn.second-candidateFirstOn.second;
        candidateMatch.Move(32,32);
        transformedCatalyst.Move(32,32);
        candidateMatch.Move(xDiff, yDiff);
        if(transformedCatalyst.Contains(candidateMatch) && candidateMatch.Contains(transformedCatalyst)){
          // further check that required cells match, too
          LifeState transformedRequired = requiredParts[catIndex];
          transformedRequired.Transform(Rotate180OddBoth);
          LifeState candidateMatchRequired = requiredParts[i];
          candidateMatchRequired.Move(xDiff, yDiff);
          if(candidateMatchRequired.Contains(transformedRequired) &&
                transformedRequired.Contains(candidateMatchRequired)){
            std::array<int,3> toPushBack({i,xDiff, yDiff});
            rotatedCatalystMatches.push_back(toPushBack);
            break;
          } else {
            std::cout << "catalysts matched up but their required parts didn't." << std::endl;
          }
        }
        if ( i + 1 == std::min(size_t(defaultCatIndex+8), catalysts.size())){
          std::cout << "Check catalyst symmetry character for catalyst number " << unorientedCatIndex;
          std::cout << " , rle " << catInput[unorientedCatIndex].rle << std::endl;
          std::cout << "After applying the transformations indicated by the symmetry character, ";
          std::cout << "one of the resulting catalysts didn't have a match when rotated by 180.";
          std::cout << "(Problematic one was number " << catIndex << ".)" << std::endl;
          exit(0);
        }
      }
    }
  }

  void Init(const char *inputFile, int nthreads) {
    begin = clock();

    std::vector<CatalystInput> inputcats;
    ReadParams(inputFile, inputcats, params);
    GenerateStates(inputcats, catalysts, requiredParts, catalystLocus, forbiddenTargets, maxMissing, transparent);
    ComputeCatalystTransformationData(inputcats);
    
    activeRegion = LifeState::Parse(params.activeRegion.c_str(), params.xActiveRegion, params.yActiveRegion);
    GenerateOffsets();

    searchArea = LifeState::SolidRect(params.searchArea[0], params.searchArea[1],
                                        params.searchArea[2], params.searchArea[3]);
    
    maxNumCatalysts = params.numCatalystsPrePost.first+params.numCatalystsPrePost.second;
    
    categoryContainer = new CategoryContainer(params.maxGen);
    fullCategoryContainer = new CategoryContainer(params.maxGen);

    for (auto & state : catalysts)
      targets.push_back(LifeTarget(state));
    
    std::cout << "last gen is " << params.lastGen << std::endl;

    // object reoccurs
    if(params.objReoccurs != std::make_pair(0,0)){
      params.filterGroups.push_back(C1);
      params.targetFilter.push_back(params.activeRegion);
      params.filterdx.push_back(params.xActiveRegion);
      params.filterdy.push_back(params.yActiveRegion);
      params.filterGen.push_back(-1);
      params.filterGenRange.push_back(params.objReoccurs);
    }

    // set up filter lists
    for (int i = 0; i < params.targetFilter.size(); i++){
      targetFilterLists.push_back(std::vector<LifeTarget>({}));
      for ( AffineTransform trans : SymmetryGroupFromEnum(params.filterGroups[i])){
        targetFilterLists[i].push_back(LifeTarget::Parse(params.targetFilter[i].c_str(),
                                        params.filterdx[i], params.filterdy[i], trans));
        // debugging purposes
        // targetFilterLists[i][targetFilterLists[i].size()-1].wanted.Print();
      }
    }

    catalystCollisionMasks = std::vector<std::vector<LifeState>>(
        catalysts.size(), std::vector<LifeState>(catalysts.size()));
    catalystReactionMasks = std::vector<LifeState>(catalysts.size());
    catalystLocusReactionMasks = std::vector<LifeState>(catalysts.size());
    catalystAvoidMasks = std::vector<LifeState>(catalysts.size());

    for (int s = 0; s < catalysts.size(); s++) {
      catalystReactionMasks[s] = catalysts[s].BigZOI();
      catalystReactionMasks[s].Transform(Rotate180OddBoth);

      LifeState nonLocus = catalysts[s];
      nonLocus.Copy(catalystLocus[s], ANDNOT);

      catalystAvoidMasks[s] = nonLocus.BigZOI();
      catalystAvoidMasks[s].Transform(Rotate180OddBoth);

      catalystLocusReactionMasks[s] = catalystLocus[s].BigZOI();
      catalystLocusReactionMasks[s].Transform(Rotate180OddBoth);
      catalystLocusReactionMasks[s].Copy(catalystAvoidMasks[s], ANDNOT);

      for (int t = 0; t < catalysts.size(); t++) {
        catalystCollisionMasks[s][t] = LoadCollisionMask(catalysts[s], catalysts[t]);
      }
    }



    current = clock();
    idx = 0;
    found = 0;
    total = 1;
    counter = 0;

    fullfound = 0;

    //int fact = 1;

    /*for (int i = 0; i < maxNumCatalysts; i++) {
      total *= params.searchArea[2];
      total *= params.searchArea[3];
      total *= catalysts.size();
      fact *= (i + 1);
    }*/

    // total /= fact;

    // std::cout << "Approximated Total: " << total << std::endl;
    // total = total / 1000000;

    //if (total == 0)
    //  total++;

    hasFilter = !params.targetFilter.empty();
    reportAll = params.fullReportFile.length() != 0;

    filterMaxGen = FilterMaxGen();

    lastReport = time(NULL);
  }

  int FilterMaxGen() {
    int maxGen = -1;

    for (int j = 0; j < targetFilterLists.size(); j++) {
      if (params.filterGen[j] > maxGen)
        maxGen = params.filterGen[j];

      if (params.filterGenRange[j].second > maxGen)
        maxGen = params.filterGenRange[j].second;
    }

    return maxGen;
  }

  void Report(const std::string& suffix) {
    std::string temp = params.outputFile;
    params.outputFile = params.outputFile + suffix + std::string(".rle");
    Report();
    params.outputFile = temp;
  }

  void Report(bool saveFile = true) const {
    float percent = ((float)idx / 10000) / (total * 1.0);
    int sec = (clock() - begin) / CLOCKS_PER_SEC + 1;
    int estimation = 0;
    int checkPerSecond = idx / (sec * 1000);

    if (percent > 0)
      estimation = (sec * 100) / percent;

    std::cout << std::setprecision(1) << std::fixed << percent << "%,"
              << idx / 1000000 << "M/" << total
              << "M, cats/total: " << categoryContainer->categories.size() << "/"
              << found;
    if (params.fullReportFile.length() != 0) {
      std::cout << ", unfiltered: " << fullCategoryContainer->categories.size() << "/"
                << fullfound;
    }
    std::cout << ", now: ";
    PrintTime(sec);
    std::cout << ", est: ";
    PrintTime(estimation);
    std::cout << ", " << std::setprecision(1) << std::fixed << checkPerSecond
              << "K/sec" << std::endl;

    // categoryContainer->Sort();
    // fullCategoryContainer->Sort();

    if (saveFile) {
      std::cout << "Saving " << params.outputFile << "... " << std::flush;

      std::ofstream catResultsFile(params.outputFile.c_str());
      catResultsFile << "x = 0, y = 0, rule = B3/S23\n";
      catResultsFile << categoryContainer->CategoriesRLE(params.maxCatSize);
      catResultsFile.close();
      std::cout << "Done!" << std::endl;

      if (params.fullReportFile.length() != 0) {
        std::cout << "Saving " << params.fullReportFile << "... " << std::flush;

        std::ofstream fullCatResultsFile(params.fullReportFile.c_str());
        fullCatResultsFile << "x = 0, y = 0, rule = B3/S23\n";
        fullCatResultsFile << fullCategoryContainer->CategoriesRLE(params.maxCatSize);
        fullCatResultsFile.close();
        std::cout << "Done!" << std::endl;
      }
    }
  }

  static void PrintTime(int sec) {
    int hr = sec / 3600;
    int min = (sec / 60) - hr * 60;
    int secs = sec - 3600 * hr - 60 * min;
    std::cout << std::setfill('0');
    std::cout << hr << ":" << std::setw(2) << min << ":" << std::setw(2) << secs;
 }

  // TODO eventually.
  bool HasForbidden(Configuration &conf, int curIter, std::pair<int,int> offset) {
    LifeState workspace;
    workspace.Join(activeRegion);
    workspace.Transform(Rotate180OddBoth);
    workspace.Move(offset.first, offset.second);
    workspace.Join(activeRegion);

    workspace.Join(conf.catalystsState); // not sure that this will quite be right

    for (int i = 0; i <= curIter + 1; i++) {
      for (int j = 0; j < maxNumCatalysts; j++) {
        for (int k = 0; k < forbiddenTargets[conf.curs[j]].size(); k++) {
          if (workspace.Contains(forbiddenTargets[conf.curs[j]][k],
                       conf.curx[j], conf.cury[j]) == true)
            return true;
        }
      }
      workspace.Step();
    }

    return false;
  }
  // TODO eventually.
  bool FilterForCurrentGenFail(LifeState &workspace) {
    for (int j = 0; j < targetFilterLists.size(); j++) {
      if (workspace.gen == params.filterGen[j] ){
        bool anyMet = false;
        for( auto filter : targetFilterLists[j]){
          anyMet = anyMet | workspace.Contains(filter);
        }
        if(!anyMet) return true;
      }
    }

    return false;
  }

  bool ValidateFilters(Configuration &conf, int catsDestroyedGen, std::pair<int,int> offset) {
    LifeState workspace = activeRegion;
    workspace.Transform(Rotate180OddBoth);
    workspace.Move(offset.first, offset.second);
    workspace.Join(activeRegion);
    workspace.Join(conf.catalystsState);

    std::vector<bool> rangeValid(params.filterGen.size(), false);

    // stop early if catalysts are destroyed
    int stopAt = params.maxGen;
    if (params.stopAfterCatsDestroyed > 0 ){
      stopAt = catsDestroyedGen + params.stopAfterCatsDestroyed;
    }
    
    for (int k = 0; k < params.filterGen.size(); k++)
      if (params.filterGen[k] >= 0)
        rangeValid[k] = true;
        // we return false early below if a single-generation filter
        // isn't met, so we don't need to keep track of those
        // via our vector rangeValid.

    for (int j = 0; j <= std::min(filterMaxGen, stopAt); j++) {
      for (int k = 0; k < params.filterGen.size(); k++) {
        if (workspace.gen == params.filterGen[k]){
          bool anyMet = false;
          for ( auto filter : targetFilterLists[k] ){
            anyMet = anyMet |  workspace.Contains(filter);
          }
          if(!anyMet)
            return false;
        }
        if (params.filterGen[k] == -1 &&
                params.filterGenRange[k].first <= workspace.gen &&
                params.filterGenRange[k].second >= workspace.gen){
          bool prevMet = rangeValid[k];
          for(LifeTarget transformedTarget :targetFilterLists[k]){
            rangeValid[k] = rangeValid[k] | workspace.Contains(transformedTarget);
          }
          if ( ! prevMet && rangeValid[k] && params.reportMatches){
            bool alreadyPresent = false;
            for (int oldMatch : matches){
              if(oldMatch == j)
                alreadyPresent = true;
            }
            if (!alreadyPresent){
              matches.push_back(j);
            }
          }
        }
      }

      workspace.Step();
    }

    for (int k = 0; k < params.filterGen.size(); k++)
      if (!rangeValid[k])
        return false;

    return true;
  }

  void ReportSolution(Configuration &conf, int successtime, const std::vector<LifeTarget> & shiftedTargets,
                      std::pair<int,int> offset){
    LifeState init;
    LifeState afterCatalyst;
    LifeState catalysts;

    LifeState workspace;
    // if reportAll - ignore filters and update fullReport
    if (reportAll) {
      workspace.Join(activeRegion);
      workspace.Transform(Rotate180OddBoth);
      workspace.Move(offset.first, offset.second);
      workspace.Join(activeRegion);
      workspace.Join(conf.catalystsState);
      init.Copy(workspace);

      workspace.Step(successtime - params.stableInterval + 2);
      afterCatalyst.Copy(workspace); // after catalyst here is post-activation but pre-recovery?

      fullfound++;

      fullCategoryContainer->Add(init, afterCatalyst, conf.catalystsState, conf,
                                 successtime - params.stableInterval + 2, 0);
      // arguments are: (LifeState &init, LifeState &afterCatalyst, LifeState &catalysts,
      //     const Configuration &conf, int firstGenSurvive,
      //     int genSurvive)
    }

    // calculate generation at which catalysts are destroyed
    int catsDestroyedGen = params.maxGen;
    if (params.stopAfterCatsDestroyed > 0){

      LifeState checkCatsDestroyed;
      checkCatsDestroyed.Join(activeRegion);
      checkCatsDestroyed.Transform(Rotate180OddBoth);
      checkCatsDestroyed.Move(offset.first, offset.second);
      checkCatsDestroyed.Join(activeRegion);
      checkCatsDestroyed.Join(conf.catalystsState);
      checkCatsDestroyed.Step(successtime);
      int absence = 0;
      
      while(checkCatsDestroyed.gen < params.maxGen+10 && absence < 10){ // 10 here more or less arbitrary.
        bool allPresent = true;
        for (auto target : shiftedTargets){
          if (!checkCatsDestroyed.Contains(target))
            allPresent = false;
        }
        if(!allPresent){
          ++absence;
        } else {
          absence = 0;
        }
        checkCatsDestroyed.Step();
      }
      catsDestroyedGen = checkCatsDestroyed.gen-absence;
      
    }

    // If has filter validate them;
    if (hasFilter) {
      if (!ValidateFilters(conf, catsDestroyedGen, offset))
        return;
    }
    // end modifications

    // if (HasForbidden(conf, successtime + 3))
    //   return;

    // If all filters validated update results
    workspace.Clear();
    workspace.Join(activeRegion);
    workspace.Transform(Rotate180OddBoth);
    workspace.Move(offset.first, offset.second);
    workspace.Join(activeRegion);
    workspace.Join(conf.catalystsState);
    init.Copy(workspace);

    workspace.Step(successtime - params.stableInterval + 2);
    afterCatalyst.Copy(workspace);

    categoryContainer->Add(init, afterCatalyst, conf.catalystsState, conf,
                           successtime - params.stableInterval + 2, 0);
    found++;
  }

  void Search() {
    Configuration config;
    config.count = 0;
    config.prePostCount = std::make_pair(0,0);
    config.transparentCount = 0;
    config.state.Join(activeRegion);
    LifeState history = config.state;
    config.postSymmetry = false;
    config.loneOffset = std::make_pair(64, 64);
    //LifeState bounds =
    //    LifeState::SolidRect(params.searchArea[0], params.searchArea[1],
    //                         params.searchArea[2], params.searchArea[3]);

    std::vector<LifeState> masks(catalysts.size());
    for (int s = 0; s < catalysts.size(); s++) {
      masks[s] = activeRegion.Convolve(catalystReactionMasks[s]);
      // masks[s].Copy(bounds, ORNOT);
    }

    LifeState required = LifeState();

    std::vector<LifeTarget> shiftedTargets(maxNumCatalysts);

    RecursiveSearch(config, history, offsets, required, masks, shiftedTargets,
                    std::array<int, MAX_CATALYSTS>(), std::array<int, MAX_CATALYSTS>(),
                    std::array<bool, MAX_CATALYSTS>(), std::array<bool, MAX_CATALYSTS>());
  }

  void
  RecursiveSearch(Configuration config, LifeState history, 
                  LifeState curOffsets,
                  const LifeState required,
                  std::vector<LifeState> masks,
                  std::vector<LifeTarget> &shiftedTargets, // This can be shared

                  std::array<int, MAX_CATALYSTS> missingTime,
                  std::array<int, MAX_CATALYSTS> recoveredTime,
                  std::array<bool, MAX_CATALYSTS> hasReacted,
                  std::array<bool, MAX_CATALYSTS> hasRecovered) {
    for (int g = config.state.gen; g < params.maxGen; g++) {
      if (config.count == 0 && g > params.lastGen){
        return;
      }
      
      if (!config.postSymmetry && g > params.lastSymInteraction){
        return;
      }

      // alternative: config.count == 0
      if (config.prePostCount.first == 0){
        if(config.postSymmetry){
          std::cout << "      Collision at gen " << g << std::endl;
        } else {
          std::cout << "Collision at gen " << g << std::endl;
        }
      }
        

      if (!config.state.Contains(required))
        return;

      for (int i = 0; i < config.count; i++) {
        // if (hasRecovered[i]) {
        //   continue;
        // }

        if (config.state.Contains(shiftedTargets[i])) {
          missingTime[i] = 0;
          recoveredTime[i] += 1;
        } else {
          hasReacted[i] = true;
          missingTime[i] += 1;
          recoveredTime[i] = 0;
        }
        if (hasReacted[i] && recoveredTime[i] > params.stableInterval)
          hasRecovered[i] = true;

        if (missingTime[i] > maxMissing[config.curs[i]])
          return;
      }

      if ( !config.postSymmetry ){
        // figure out which offsets will interact next gen
        LifeState offsetsThatInteract = curOffsets;
        
        // eliminate what you can based off of bounding boxes [worth doing?]
        std::array<int,4> xyBounds = config.state.XYBounds();
        LifeState box = LifeState::SolidRect(2*xyBounds[0]-2,2*xyBounds[1]-2, 2*xyBounds[2]-2*xyBounds[0]+4,
                    2*xyBounds[3]-2*xyBounds[1]+4);
        offsetsThatInteract.Copy(box, AND);

        
        LifeState neighborhood = config.state.Convolve(LifeState::SolidRect(-1,-1,3,3));
        offsetsThatInteract.Copy(neighborhood.Convolve(neighborhood), AND);

        // remove those offsets that interact
        curOffsets.Copy(offsetsThatInteract, ANDNOT);
        
        while( !offsetsThatInteract.IsEmpty()){
          std::pair<int,int> offset = offsetsThatInteract.FirstOn();
          offsetsThatInteract.Erase(offset.first, offset.second);

          // could do a brief lookahead to see if required cells will die soon.
          if(config.prePostCount.first > 0 ){
            LifeState lookahead = config.state;
            lookahead.Transform(Rotate180OddBoth);
            lookahead.Move(offset.first, offset.second);
            lookahead.Join(config.state);
            lookahead.Step();
            lookahead.Step();
            lookahead.Step();
            lookahead.Step();
            if(!lookahead.Contains(required)){
              continue; // is this continue a problem?
              // it shouldn't be, because we erase first thing.
            }
          }

          Configuration newConfig = config;
          newConfig.postSymmetry = true;
          newConfig.loneOffset = offset;
          newConfig.state.Transform(Rotate180OddBoth);
          newConfig.state.Move(offset.first, offset.second);
          newConfig.state.Join(config.state);

          newConfig.catalystsState.Transform(Rotate180OddBoth);
          newConfig.catalystsState.Move(offset.first, offset.second);
          newConfig.catalystsState.Join(config.catalystsState);

          newConfig.catsPreSymmetry.Transform(Rotate180OddBoth);
          newConfig.catsPreSymmetry.Move(offset.first, offset.second);
          newConfig.catsPreSymmetry.Join(config.catsPreSymmetry);

          LifeState newHistory = history;
          newHistory.Transform(Rotate180OddBoth);
          newHistory.Move(offset.first, offset.second);
          newHistory.Join(history);

          std::vector<LifeState> newMasks = masks;
          for (int i = 0; i < masks.size(); ++i ){
            // if you're checking masks[s] at (x,y), you should check
            // masks[rotatedS] at (rotatedX0-x,rotatedY0-y) too.
            // ie masks[s] should be combined with {masks[rotatedS],
            // shifted by -rotatedX0, -rotatedY0, rotated by 180 }
            int rotatedS = rotatedCatalystMatches[i][0];
            int rotatedX0 = rotatedCatalystMatches[i][1];
            int rotatedY0 = rotatedCatalystMatches[i][2];
            LifeState transformedMask = masks[rotatedS];
            transformedMask.Transform(-rotatedX0, -rotatedY0, Rotate180OddBoth);
            // prohibited by either => prohibited, so this really is OR not AND
            newMasks[i].Copy(transformedMask, OR); 
          }
          
          LifeState newOffsets = LifeState::SolidRect(offset.first, offset.second, 1,1);

          // don't need to change targets or required.
          RecursiveSearch(newConfig, newHistory, newOffsets, required, newMasks,
                    shiftedTargets, missingTime, recoveredTime, hasReacted,
                    hasRecovered);
        
        } // end of going thru offsets that interact.

      } // end of checking which interact


      // Try adding a catalyst
      bool donePlacing = config.postSymmetry ?  config.prePostCount.second == params.numCatalystsPrePost.second :
                                                    config.prePostCount.first == params.numCatalystsPrePost.first;
      
      if (config.state.gen >= params.startGen && !donePlacing) {
        LifeState newcells = config.state;
        newcells.Copy(history, ANDNOT);

        if (!newcells.IsEmpty()) {
          // for (int s = 0; s < catalysts.size(); s++) {
          //   LifeState hitLocations = newcells.Convolve(catalystAvoidMasks[s]);
          //   masks[s].Join(hitLocations);
          // }

          for (int s = 0; s < catalysts.size(); s++) {
            
            if (transparent[s] && config.transparentCount == params.numTransparent)
              continue;

            LifeState newPlacements = catalystReactionMasks[s].Convolve(newcells);
            newPlacements.Copy(masks[s], ANDNOT);
            newPlacements.Copy(searchArea, AND);

            while (!newPlacements.IsEmpty()) {
              // Do the placement
              auto newPlacement = newPlacements.FirstOn();

              LifeState offsetCatValid = curOffsets;
              if( !config.postSymmetry ){
                // rotatedCatalystMatches[s] := rotatedS, rotatedX, rotatedY
                // want all v's such that placing catalysts[rotatedS] @ rotVec-newPlaceVec+v is ok.
                int rotatedS = rotatedCatalystMatches[s][0];
                int rotatedX0 = rotatedCatalystMatches[s][1];
                int rotatedY0 = rotatedCatalystMatches[s][2];
                LifeState rotatedCatMask = masks[rotatedS];
                rotatedCatMask.Move(newPlacement.first-rotatedX0, newPlacement.second-rotatedY0);
                offsetCatValid.Copy(rotatedCatMask, ANDNOT);
                if (offsetCatValid.IsEmpty()){
                  masks[s].Set(newPlacement.first, newPlacement.second);
                  newPlacements.Erase(newPlacement.first, newPlacement.second);
                  continue;
                }
              }

              if (config.count == 0 && !params.quietMode) {
                if (config.postSymmetry)
                  std::cout << "     ";
                std::cout << "Placing catalyst " << s << " at " \
                          << newPlacement.first << ", " << newPlacement.second \
                          << std::endl;
                // newPlacements.Print();
              }

              Configuration newConfig = config;
              newConfig.count += 1;
              if( config.postSymmetry){
                newConfig.prePostCount.second += 1;
              } else {
                newConfig.prePostCount.first += 1;
              }
              newConfig.curx[config.count] = newPlacement.first;
              newConfig.cury[config.count] = newPlacement.second;
              newConfig.curs[config.count] = s;
              if (transparent[s])
                newConfig.transparentCount++;

              LifeState shiftedCatalyst = catalysts[s];
              shiftedCatalyst.Move(newPlacement.first, newPlacement.second);
              LifeState newHistory = history;
              if(!config.postSymmetry){
                newConfig.catalystsState.Join(shiftedCatalyst);
                newConfig.state.Join(shiftedCatalyst);
                newHistory.Join(shiftedCatalyst);

                newConfig.catsPreSymmetry.Join(shiftedCatalyst);
              } else {
                assert(config.loneOffset == curOffsets.FirstOn());
                LifeState symCatalyst = shiftedCatalyst;
                symCatalyst.Transform(Rotate180OddBoth);
                symCatalyst.Move(config.loneOffset.first, config.loneOffset.second);
                symCatalyst.Join(shiftedCatalyst);
                newConfig.catalystsState.Join(symCatalyst);
                newConfig.state.Join(symCatalyst);
                newHistory.Join(symCatalyst);

                newConfig.catsPostSymmetry.Join(symCatalyst);
              }

              LifeState newRequired = required;
              newRequired.Join(requiredParts[s], newPlacement.first, newPlacement.second);

              if (newConfig.prePostCount.second != params.numCatalystsPrePost.second) {
                LifeState lookahead = newConfig.state;
                lookahead.Step();
                lookahead.Step();
                lookahead.Step();
                lookahead.Step();
                if (!lookahead.Contains(newRequired)) {
                  masks[s].Set(newPlacement.first, newPlacement.second);
                  newPlacements.Erase(newPlacement.first, newPlacement.second);
                  continue;
                }
              }

              std::vector<LifeState> newMasks = masks;

              shiftedTargets[config.count] = LifeTarget(shiftedCatalyst);

              LifeState bounds;
              if (params.maxW != -1) {
                bounds = LifeState::SolidRect(newPlacement.first - params.maxW,
                                              newPlacement.second - params.maxH,
                                              2 * params.maxW - 1,
                                              2 * params.maxH - 1);
              }

              // If we just placed the last catalyst, don't bother
              if (newConfig.prePostCount.second != params.numCatalystsPrePost.second) {
                for (int t = 0; t < catalysts.size(); t++) {
                  newMasks[t].Join(catalystCollisionMasks[s][t],
                                   newPlacement.first, newPlacement.second);

                  if (params.maxW != -1) {
                    newMasks[t].Copy(bounds, ORNOT);
                  }
                }
              }

              RecursiveSearch(newConfig, newHistory, offsetCatValid, newRequired, newMasks, shiftedTargets, missingTime,
                              recoveredTime, hasReacted, hasRecovered);
              
              /*   RecursiveSearch(Configuration config, LifeState history, 
                  LifeState curOffsets,
                  const LifeState required,
                  std::vector<LifeState> masks,
                  std::vector<LifeTarget> &shiftedTargets, // This can be shared

                  std::array<int, MAX_CATALYSTS> missingTime,
                  std::array<int, MAX_CATALYSTS> recoveredTime,
                  std::array<bool, MAX_CATALYSTS> hasReacted,
                  std::array<bool, MAX_CATALYSTS> hasRecovered)*/

              masks[s].Set(newPlacement.first, newPlacement.second); 
              newPlacements.Erase(newPlacement.first, newPlacement.second); // is there a reason why this
              // is after the recursive call, not before?

            } // end while placements not empty.
          } // end for catalysts
        } // end newcells
      } // end of "if there's more catalysts to place"

      // Still block the locations that are hit too early
      if (config.state.gen < params.startGen) {
        for (int s = 0; s < catalysts.size(); s++) {
          LifeState hitLocations = config.state.Convolve(catalystReactionMasks[s]);
          masks[s].Join(hitLocations);
        }
      }
      if (config.prePostCount.second == params.numCatalystsPrePost.second) {
        bool allRecovered = true;
        for (int i = 0; i < config.count; i++) {
          if (!hasRecovered[i] || missingTime[i] > 0) {
            allRecovered = false;
          }
        }
        if(allRecovered) {
          ReportSolution(config, g, shiftedTargets, config.loneOffset);
          return;
        }
      }

      // save at least every 15 minutes or so
      if ( (config.count == 0 && difftime(time(NULL), lastReport) > 10) || \
              (config.count == 1 && difftime(time(NULL), lastReport) > 900) ){
        Report();
        lastReport = time(NULL);
      }
      history.Copy(config.state, OR);
      config.state.Step();
    }
  }
};

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cout << "Usage CatForce.exe <in file> <nthreads>" << std::endl;
    exit(0);
  }
  int nthreads = 1;
  if (argc > 2) {
    nthreads = atoi(argv[2]);
  }
  // omp_set_num_threads(nthreads);
  std::cout << "Input: " << argv[1] << std::endl
            << "Initializing please wait..." << std::endl;

  CatalystSearcher searcher;
  searcher.Init(argv[1], nthreads);

  clock_t initialized = clock();
  printf("Total elapsed CPU time (not wallclock if nthreads>1): %f seconds\n",
         (double) (initialized - searcher.begin) / CLOCKS_PER_SEC);
  
  std::cout << std::endl
            << "Initialization finished, searching..." << std::endl
            << std::endl;

  searcher.Search();
  // Print report one final time (update files with the final results).
  searcher.Report();
  printf("\n\nFINISH\n");
  clock_t end = clock();
  printf("Total elapsed CPU time (not wallclock if nthreads>1): %f seconds\n",
         (double)(end - searcher.begin) / CLOCKS_PER_SEC);
  if(searcher.matches.size() != 0){
    std::cout << "Ranged filters matched at generations ";
    for(int i = 0 ; i < searcher.matches.size(); ++i){
      std::cout << searcher.matches[i];
      if (i + 1 < searcher.matches.size()){
        std::cout << " ";
      }
    }
    std::cout << std::endl;
  }
}