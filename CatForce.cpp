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
  unsigned maxGen;
  unsigned numCatalysts;
  unsigned numTransparent;
  unsigned stableInterval;
  std::string pat;
  int xPat;
  int yPat;
  unsigned startGen;
  unsigned lastGen;
  std::string outputFile;
  std::string fullReportFile;
  int searchArea[4]{};
  int maxW;
  int maxH;
  StaticSymmetry symmetryEnum;

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
  bool firstTransparent;

  std::vector<std::pair<int, int>> filterGenRange;


  std::tuple<std::string, int, int> alsoRequired;

  int maxCatSize;

  SearchParams() {
    maxGen = 250;
    numCatalysts = 2;
    stableInterval = 15;
    pat = "";
    searchArea[0] = -10;
    searchArea[1] = 0;
    searchArea[2] = 20;
    searchArea[3] = 20;
    xPat = 0;
    yPat = 0;
    startGen = 1;
    lastGen = 100;
    outputFile = "results.rle";
    maxW = -1;
    maxH = -1;
    symmetryEnum = StaticSymmetry::C1;
    maxCatSize = -1;
    fullReportFile = "";

    alsoRequired = std::make_tuple("", 0,0);
    
    filterGroups = {};
    stopAfterCatsDestroyed = -1;

    quietMode = false;
    reportMatches = false;
    firstTransparent = false;
  }
};

class CatalystInput {
public:
  std::string rle;
  unsigned maxDisappear;
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
  bool mustInclude;

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
    mustInclude = false;

    unsigned argi = 6;

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
      } else if (elems[argi] == "mustinclude") {
        mustInclude = true;
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

  std::string numCat = "num-catalyst";
  std::string numTransp = "num-transparent";
  std::string stable = "stable-interval";
  std::string area = "search-area";
  std::string pat = "pat";
  std::string outputFile = "output";
  std::string filter = "filter";
  std::string maxWH = "fit-in-width-height";
  std::string maxCatSize = "max-category-size";
  std::string fullReport = "full-report";

  std::string symmetry = "symmetry";
  std::string alsoRequired = "also-required";

  std::string stopAfterCatsDestroyed = "stop-after-cats-destroyed";
  std::string reportMatches = "report-matches";
  std::string quiet = "quiet-mode";
  std::string firstTranps = "first-transparent";

  std::string line;
  bool hasLastGen = false;

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

      if (elems[0] == numCat)
        params.numCatalysts = atoi(elems[1].c_str());

      if (elems[0] == numTransp)
        params.numTransparent = atoi(elems[1].c_str());

      if (elems[0] == stable)
        params.stableInterval = atoi(elems[1].c_str());

      if (elems[0] == pat) {
        params.pat = elems[1];

        if (elems.size() > 3) {
          params.xPat = atoi(elems[2].c_str());
          params.yPat = atoi(elems[3].c_str());
        }
      }

      if (elems[0] == area) {
        params.searchArea[0] = atoi(elems[1].c_str());
        params.searchArea[1] = atoi(elems[2].c_str());
        params.searchArea[2] = atoi(elems[3].c_str());
        params.searchArea[3] = atoi(elems[4].c_str());
      }

      if (elems[0] == startGen)
        params.startGen = atoi(elems[1].c_str());

      if (elems[0] == lastGen) {
        params.lastGen = atoi(elems[1].c_str());
        hasLastGen = true;
      }

      if (elems[0] == outputFile) {
        params.outputFile = elems[1];

        for (unsigned i = 2; i < elems.size(); i++) {
          params.outputFile.append(" ");
          params.outputFile.append(elems[i]);
        }
      }

      if (elems[0] == fullReport) {
        params.fullReportFile = elems[1];

        for (unsigned i = 2; i < elems.size(); i++) {
          params.fullReportFile.append(" ");
          params.fullReportFile.append(elems[i]);
        }
      }

      if (elems[0] == filter) {
        std::vector<std::string> rangeElems;
        split(elems[1], '-', rangeElems);

        if (rangeElems.size() == 1) {
          params.filterGen.push_back(atoi(elems[1].c_str()));
          params.filterGenRange.emplace_back(-1, -1);
        } else {
          unsigned minGen = atoi(rangeElems[0].c_str());
          unsigned maxGen = atoi(rangeElems[1].c_str());

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

      if (elems[0] == alsoRequired){
        params.alsoRequired = std::make_tuple(elems[1], atoi(elems[2].c_str()), atoi(elems[3].c_str()));
      }

      if (elems[0] == symmetry) {
        if ( elems[1] == "C1" || SymmetryEnumFromString(elems[1]) != StaticSymmetry::C1){
          params.symmetryEnum = SymmetryEnumFromString(elems[1]);
        } else {
          std::cout << "Couldn't parse symmetry option " << elems[1] << std::endl;
          exit(0);
        }
      }

      if (elems[0] == stopAfterCatsDestroyed){
        params.stopAfterCatsDestroyed = atoi(elems[1].c_str());
      }

      if (elems[0] == reportMatches){
        params.reportMatches = true;
      }


      if (elems[0] == firstTranps){
        params.firstTransparent = true;
      }


      if (elems[0] == quiet){
        params.quietMode = true;
      }

    } catch (const std::exception &ex) {
    }
  }
  if(!hasLastGen)
    params.lastGen = params.maxGen - 1;

  if (params.pat.length() == 0) {
    std::cout << "Did not read any pattern!" << std::endl;
    exit(1);
  }
  if (catalysts.empty()) {
    std::cout << "Did not read any catalysts!" << std::endl;
    exit(1);
  }
}

class CatalystData {
public:
  LifeState state;
  // LifeTarget target;
  LifeState reactionMask;
  unsigned maxDisappear;
  std::vector<LifeTarget> forbidden;
  LifeState required;
  LifeState locus;
  bool transparent;
  bool mustInclude;

  static std::vector<CatalystData> FromInput(CatalystInput &input);
};

std::vector<CatalystData> CatalystData::FromInput(CatalystInput &input) {
  std::vector<AffineTransform> trans;
  CharToTransVec(input.symmType, trans);

  const char *rle = input.rle.c_str();

  std::vector<CatalystData> results;

  for (auto &tran : trans) {
    LifeState pat = LifeState::Parse(rle, input.centerX, input.centerY, tran);

    CatalystData result;

    result.state = pat;
    result.reactionMask = pat.BigZOI();
    result.reactionMask.Transform(Rotate180OddBoth);

    // result.target = LifeTarget(pat);
    result.maxDisappear = input.maxDisappear;

    for (unsigned k = 0; k < input.forbiddenRLE.size(); k++) {
      result.forbidden.push_back(LifeTarget::Parse(input.forbiddenRLE[k].c_str(),
                                                   input.forbiddenXY[k].first,
                                                   input.forbiddenXY[k].second, tran));
    }

    if (input.requiredRLE != "") {
      result.required = LifeState::Parse(input.requiredRLE.c_str(),
                                         input.requiredXY.first,
                                         input.requiredXY.second, tran);
    }

    if (input.locusRLE != "") {
      result.locus = LifeState::Parse(input.locusRLE.c_str(),
                                      input.locusXY.first,
                                      input.locusXY.second, tran);
    } else {
      result.locus = pat;
    }

    result.transparent = input.transparent;
    result.mustInclude = input.mustInclude;

    results.push_back(result);
  }
  return results;
}

struct Configuration {
  unsigned count;
  unsigned transparentCount;
  unsigned mustIncludeCount;
  std::array<int, MAX_CATALYSTS> curx;
  std::array<int, MAX_CATALYSTS> cury;
  std::array<int, MAX_CATALYSTS> curs;
  // int minIter;
  LifeState state;
  LifeState catalystsState;
  int firstReaction;
  bool anyTransparent;

  LifeState startingCatalysts;
  // std::vector<LifeTarget> shiftedTargets;
};

// Fix a, what positions of b causes a collision?
LifeState CollisionMask(const LifeState &a, const LifeState &b) {
  unsigned popsum = a.GetPop() + b.GetPop();

  LifeState mask;
  for (unsigned x = 0; x < N; x++) {
    for (unsigned y = 0; y < 64; y++) {
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

LifeState LoadCollisionMask(const CatalystData &a, const CatalystData &b) {
  std::stringstream ss;
  ss << "masks/mask-" << a.state.GetHash() << "-" << b.state.GetHash();
  std::string fname = ss.str();

  std::ifstream infile;
  infile.open(fname.c_str(), std::ifstream::in);
  if (!infile.good()) {
    LifeState mask = CollisionMask(a.state, b.state);
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

std::string GetRLE(const std::vector<std::vector<bool>> &life2d) {
  if (life2d.empty())
    return "";

  if (life2d[0].empty())
    return "";

  std::stringstream result;

  unsigned eol_count = 0;

  for (unsigned j = 0; j < life2d[0].size(); j++) {
    bool last_val = life2d[0][j];
    unsigned run_count = 0;

    for (const auto & i : life2d) {
      bool val = i[j];

      // Flush linefeeds if we find a live cell
      if (val && eol_count > 0) {
        if (eol_count > 1)
          result << eol_count;

        result << "$";

        eol_count = 0;
      }

      // Flush current run if val changes
      if (val == !last_val) {
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
    if (last_val) {
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
  std::vector<std::vector<bool>> vec(N, std::vector<bool>(N));

  for (unsigned j = 0; j < N; j++)
    for (unsigned i = 0; i < N; i++)
      vec[i][j] = s.GetCell(i - 32, j - 32) == 1;

  return GetRLE(vec);
}

class SearchResult {
public:
  // Saved for the report
  LifeState init;

  // iters state in form of integers
  // std::vector<int> params;
  unsigned maxGenSurvive;
  unsigned firstGenSurvive;

  SearchResult(LifeState &initState, const Configuration &conf,
               unsigned firstGenSurviveIn, unsigned genSurvive) {
    init.Copy(initState);

    maxGenSurvive = genSurvive;
    firstGenSurvive = firstGenSurviveIn;
  }

  void Print() {
    std::cout << "start:" << firstGenSurvive;
    std::cout << ", finish:" << maxGenSurvive << ", params: ";

    // for (int param : params)
    //   std::cout << param << ",";

    std::cout << std::endl;
  }
};

class Category {
private:
  unsigned catDelta;
  int maxgen;
  uint64_t hash;

public:
  LifeState categoryKey;
  std::vector<SearchResult> results;

  Category(LifeState &catalystRemoved, SearchResult &firstResult,
           unsigned catDeltaIn, unsigned maxGen) {
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

    for (unsigned i = 0; i < catDelta; i++) {
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
    const unsigned Dist = 36 + 64;

    unsigned howmany = results.size();

    if (maxCatSize != -1)
      howmany = std::min(howmany, (unsigned)maxCatSize);

    unsigned width = Dist * howmany;
    unsigned height = Dist;

    std::vector<std::vector<bool>> vec(width, std::vector<bool>(height));

    for (unsigned l = 0; l < howmany; l++)
      for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++)
          vec[Dist * l + i][j] = results[l].init.GetCell(i - 32, j - 32) == 1;

    return GetRLE(vec);
  }
};

class CategoryContainer {
public:
  std::vector<Category *> categories;
  unsigned catDelta;
  unsigned maxgen;

  explicit CategoryContainer(unsigned maxGen) {
    catDelta = 14;
    maxgen = maxGen + catDelta;
  }

  CategoryContainer(unsigned cats, unsigned maxGen) {
    catDelta = cats;
    maxgen = maxGen + catDelta;
  }

  void Add(LifeState &init, const LifeState &afterCatalyst, const LifeState &catalysts,
           const Configuration &conf, unsigned firstGenSurvive,
           unsigned genSurvive) {
    LifeState result;

    result.Copy(afterCatalyst);
    result.Copy(catalysts, XOR);

    result.Step(maxgen - result.gen);
    uint64_t hash = result.GetHash();

    result.Clear();
    result.Copy(afterCatalyst);
    result.Copy(catalysts, XOR);

    for (auto & category: categories) {
      if (category->BelongsTo(result, hash)) {
          SearchResult r(init, conf, firstGenSurvive, genSurvive);
          category->Add(r);
          return;
      }
    }

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
  LifeState pat;
  std::vector<std::vector<LifeTarget>> targetFilterLists;
  std::vector<CatalystData> catalysts;
  std::vector<std::vector<LifeState>> catalystCollisionMasks;

  clock_t current{};
  long long idx{};
  unsigned found{};
  unsigned fullfound{};
  long long total{};
  unsigned short int counter{};

  CategoryContainer *categoryContainer{};
  CategoryContainer *fullCategoryContainer{};

  bool hasFilter{};
  bool reportAll{};

  int filterMaxGen{};

  std::vector<unsigned> matches;

  time_t lastReport;
  void Init(const char *inputFile) {
    begin = clock();

    std::vector<CatalystInput> inputcats;
    ReadParams(inputFile, inputcats, params);

    for (auto &input : inputcats) {
      std::vector<CatalystData> newcats = CatalystData::FromInput(input);
      catalysts.insert(catalysts.end(), newcats.begin(), newcats.end());
    }
    bool hasMustInclude = false;
    for (auto &cat : catalysts) {
      hasMustInclude = hasMustInclude || cat.mustInclude;
    }

    if (!hasMustInclude) {
      for (auto &cat : catalysts) {
        cat.mustInclude = true;
      }
    }

    pat = LifeState::Parse(params.pat.c_str(), params.xPat, params.yPat);
    categoryContainer = new CategoryContainer(params.maxGen);
    fullCategoryContainer = new CategoryContainer(params.maxGen);

    // modification
    for (unsigned i = 0; i < params.targetFilter.size(); i++){
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

    for (unsigned s = 0; s < catalysts.size(); s++) {
      // LifeState nonLocus = catalysts[s];
      // nonLocus.Copy(catalystLocus[s], ANDNOT);

      // catalystAvoidMasks[s] = nonLocus.BigZOI();
      // catalystAvoidMasks[s].Transform(Rotate180OddBoth);

      // catalystLocusReactionMasks[s] = catalystLocus[s].BigZOI();
      // catalystLocusReactionMasks[s].Transform(Rotate180OddBoth);
      // catalystLocusReactionMasks[s].Copy(catalystAvoidMasks[s], ANDNOT);

      for (unsigned t = 0; t < catalysts.size(); t++) {
        catalystCollisionMasks[s][t] = LoadCollisionMask(catalysts[s], catalysts[t]);
      }
    }

    current = clock();
    found = 0;

    hasFilter = !params.targetFilter.empty();
    reportAll = params.fullReportFile.length() != 0;

    filterMaxGen = FilterMaxGen();

    lastReport = time(NULL);
  }

  int FilterMaxGen() {
    int maxGen = 0;

    for (unsigned j = 0; j < params.targetFilter.size(); j++) {
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
    unsigned sec = (clock() - begin) / CLOCKS_PER_SEC + 1;
    // unsigned checkPerSecond = idx / (sec * 1000);

    std::cout << "results: " << categoryContainer->categories.size() << "/"
              << found;
    if (params.fullReportFile.length() != 0) {
      std::cout << ", unfiltered: " << fullCategoryContainer->categories.size() << "/"
                << fullfound;
    }
    std::cout << ", now: ";
    PrintTime(sec);
    std::cout << std::endl;
    // std::cout << ", " << std::setprecision(1) << std::fixed << checkPerSecond
    //           << "K/sec" << std::endl;

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

  static void PrintTime(unsigned sec) {
    unsigned hr = sec / 3600;
    unsigned min = (sec / 60) - hr * 60;
    unsigned secs = sec - 3600 * hr - 60 * min;
    std::cout << std::setfill('0');
    std::cout << hr << ":" << std::setw(2) << min << ":" << std::setw(2) << secs;
 }

  bool HasForbidden(Configuration &conf, unsigned curIter) {
    LifeState workspace;
    workspace.Join(conf.startingCatalysts);
    workspace.JoinWSymChain(pat, params.symmetryEnum);

    for (unsigned i = 0; i <= curIter + 1; i++) {
      for (unsigned j = 0; j < catalysts.size(); j++) {
        for (unsigned k = 0; k < catalysts[conf.curs[j]].forbidden.size(); k++) {
          if (workspace.Contains(catalysts[conf.curs[j]].forbidden[k],
                                 conf.curx[j], conf.cury[j]))
            return true;
        }
      }
      workspace.Step();
    }

    return false;
  }

  bool FilterForCurrentGenFail(LifeState &workspace) {
    for (unsigned j = 0; j < targetFilterLists.size(); j++) {
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

  bool ValidateFilters(Configuration &conf, int catsDestroyedGen) {
    LifeState workspace;
    workspace.JoinWSymChain(pat, params.symmetryEnum);
    workspace.Join(conf.catalystsState);

    std::vector<bool> rangeValid(params.filterGen.size(), false);

    // stop early if catalysts are destroyed
    int stopAt = params.maxGen;
    if (params.stopAfterCatsDestroyed > 0 ){
      stopAt = catsDestroyedGen + params.stopAfterCatsDestroyed;
    }
    
    for (unsigned k = 0; k < params.filterGen.size(); k++)
      if (params.filterGen[k] >= 0)
        rangeValid[k] = true;
        // we return false early below if a single-generation filter
        // isn't met, so we don't need to keep track of those
        // via our vector rangeValid.

    for (unsigned j = 0; j <= std::min(filterMaxGen, stopAt); j++) {
      for (unsigned k = 0; k < params.filterGen.size(); k++) {
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
            for (unsigned oldMatch : matches){
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

    for (unsigned k = 0; k < params.filterGen.size(); k++)
      if (!rangeValid[k])
        return false;

    return true;
  }

  void ReportSolution(Configuration &conf, unsigned successtime, std::vector<LifeTarget> & shiftedTargets){
    LifeState init;
    LifeState afterCatalyst;
    LifeState catalysts;

    LifeState workspace;
    // if reportAll - ignore filters and update fullReport
    if (reportAll) {
      workspace.Join(conf.startingCatalysts);
      workspace.JoinWSymChain(pat, params.symmetryEnum);
      init.Copy(workspace);

      workspace.Step(successtime - params.stableInterval + 2);
      afterCatalyst.Copy(workspace);

      fullfound++;

      fullCategoryContainer->Add(init, afterCatalyst, conf.startingCatalysts, conf,
                                 successtime - params.stableInterval + 2, 0);
    }

    // calculate generation at which catalysts are destroyed
    int catsDestroyedGen = params.maxGen;
    if (params.stopAfterCatsDestroyed > 0){

      LifeState checkCatsDestroyed;
      checkCatsDestroyed.Join(conf.catalystsState);
      checkCatsDestroyed.JoinWSymChain(pat,  params.symmetryEnum);
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
      if (!ValidateFilters(conf, catsDestroyedGen))
        return;
    }
    // end modifications

    // if (HasForbidden(conf, successtime + 3))
    //   return;

    // If all filters validated update results
    workspace.Clear();
    workspace.JoinWSymChain(pat, params.symmetryEnum);
    workspace.Join(conf.startingCatalysts);
    init.Copy(workspace);

    workspace.Step(successtime - params.stableInterval + 2);
    afterCatalyst.Copy(workspace);

    categoryContainer->Add(init, afterCatalyst, conf.startingCatalysts, conf,
                           successtime - params.stableInterval + 2, 0);
    found++;
  }

  void Search() {
    Configuration config;
    config.firstReaction = 1000; // large number.
    config.anyTransparent = false;
    config.count = 0;
    config.transparentCount = 0;
    config.mustIncludeCount = 0;
    config.state.JoinWSymChain(pat, params.symmetryEnum);
    LifeState history = config.state;

    LifeState bounds =
        LifeState::SolidRect(params.searchArea[0], params.searchArea[1],
                             params.searchArea[2], params.searchArea[3]);

    std::vector<LifeState> masks(catalysts.size());
    for (unsigned s = 0; s < catalysts.size(); s++) {
      masks[s] = config.state.Convolve(catalysts[s].reactionMask);
      masks[s].Copy(bounds, ORNOT);
    }

    LifeState required = LifeState();
    if(std::get<0>(params.alsoRequired) != ""){

      LifeState alsoReq = LifeState::Parse(std::get<0>(params.alsoRequired).c_str(), std::get<1>(params.alsoRequired),
                                                                                        std::get<2>(params.alsoRequired));
      if (! config.state.Contains(alsoReq)){
        std::cout << "also-required parameter was passed, yet said pattern isn't in the input RLE." << std::endl;
        exit(0);
      }
      required.Join(alsoReq);
    }

    std::vector<LifeTarget> shiftedTargets(params.numCatalysts);

    RecursiveSearch(config, history, required, masks, shiftedTargets,
                    std::array<unsigned, MAX_CATALYSTS>(), std::array<unsigned, MAX_CATALYSTS>(),
                    std::array<bool, MAX_CATALYSTS>(), std::array<bool, MAX_CATALYSTS>());
  }

  void
  RecursiveSearch(Configuration config, LifeState history, const LifeState required,
                  std::vector<LifeState> masks,
                  std::vector<LifeTarget> &shiftedTargets, // This can be shared

                  std::array<unsigned, MAX_CATALYSTS> missingTime,
                  std::array<unsigned, MAX_CATALYSTS> recoveredTime,
                  std::array<bool, MAX_CATALYSTS> hasReacted,
                  std::array<bool, MAX_CATALYSTS> hasRecovered) {

    for (unsigned g = config.state.gen; g < params.maxGen; g++) {
      if (config.count == 0 && g > params.lastGen)
        return;

      if (config.count == 0) {
        std::cout << "Collision at gen " << g << std::endl;
      }

      if (g == config.firstReaction + params.stableInterval && params.firstTransparent  && !config.anyTransparent){
        return;
      }

      if (!config.state.Contains(required))
        return;
      
      // if we've placed all the catalysts and none are marked as transparent, return.
      if (params.firstTransparent && config.count == params.numCatalysts){
        bool anyTransparent = false;
        for(int i = 0; i < config.count; ++i){
          anyTransparent = anyTransparent | catalysts[i].transparent;
        }
        if (!anyTransparent){
          return;
        }
      }

      for (unsigned i = 0; i < config.count; i++) {
        // if (hasRecovered[i]) {
        //   continue;
        // }

        if(params.firstTransparent && catalysts[i].transparent){
           if(config.firstReaction < params.maxGen && g < config.firstReaction + params.stableInterval
                                                    && config.state.AreDisjoint(shiftedTargets[i].wanted) ){
            config.anyTransparent = true;
          }
        }

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

        if (missingTime[i] > catalysts[config.curs[i]].maxDisappear)
          return;
      }

      // Try adding a catalyst
      if (config.state.gen >= params.startGen && config.count != params.numCatalysts) {
        LifeState newcells = config.state;
        newcells.Copy(history, ANDNOT);

        if (!newcells.IsEmpty()) {
          // for (unsigned s = 0; s < catalysts.size(); s++) {
          //   LifeState hitLocations = newcells.Convolve(catalystAvoidMasks[s]);
          //   masks[s].Join(hitLocations);
          // }

          for (unsigned s = 0; s < catalysts.size(); s++) {
            if (config.transparentCount == params.numTransparent && catalysts[s].transparent)
              continue;
            if (config.count == params.numCatalysts - 1 && config.mustIncludeCount == 0 && !catalysts[s].mustInclude)
              continue;

            LifeState newPlacements = catalysts[s].reactionMask.Convolve(newcells);
            newPlacements.Copy(masks[s], ANDNOT);

            while (!newPlacements.IsEmpty()) {
              // Do the placement
              auto newPlacement = newPlacements.FirstOn();

              if (config.count == 0 && !params.quietMode) {
                std::cout << "Placing catalyst " << s << " at "
                          << newPlacement.first << ", " << newPlacement.second
                          << std::endl;
              }

              Configuration newConfig = config;
              newConfig.count += 1;
              newConfig.curx[config.count] = newPlacement.first;
              newConfig.cury[config.count] = newPlacement.second;
              newConfig.curs[config.count] = s;
              if (catalysts[s].transparent)
                newConfig.transparentCount++;
              if (catalysts[s].mustInclude)
                newConfig.mustIncludeCount++;

              LifeState shiftedCatalyst = catalysts[s].state;
              shiftedCatalyst.Move(newPlacement.first, newPlacement.second);
              shiftedTargets[config.count] = LifeTarget(shiftedCatalyst);

              LifeState symCatalyst;
              symCatalyst.JoinWSymChain(shiftedCatalyst, params.symmetryEnum);
              newConfig.startingCatalysts.Join(symCatalyst);
              newConfig.catalystsState.Join(symCatalyst);
              newConfig.state.Join(symCatalyst);

              LifeState newHistory = history;
              newHistory.Join(symCatalyst);

              LifeState newRequired = required;
              newRequired.Join(catalysts[s].required, newPlacement.first, newPlacement.second);

              if (newConfig.count != params.numCatalysts) {
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

              LifeState bounds;
              if (params.maxW != -1) {
                bounds = LifeState::SolidRect(newPlacement.first - params.maxW,
                                              newPlacement.second - params.maxH,
                                              2 * params.maxW - 1,
                                              2 * params.maxH - 1);
              }

              if(config.count == 0)
                newConfig.firstReaction = g;

              // If we just placed the last catalyst, don't bother
              if (newConfig.count != params.numCatalysts) {
                for (unsigned t = 0; t < catalysts.size(); t++) {
                  newMasks[t].Join(catalystCollisionMasks[s][t],
                                   newPlacement.first, newPlacement.second);

                  if (params.maxW != -1) {
                    newMasks[t].Copy(bounds, ORNOT);
                  }
                }
              }

              RecursiveSearch(newConfig, newHistory, newRequired, newMasks, shiftedTargets, missingTime,
                              recoveredTime, hasReacted, hasRecovered);

              masks[s].Set(newPlacement.first, newPlacement.second);
              newPlacements.Erase(newPlacement.first, newPlacement.second);
            }
          }
        }
      }

      // Still block the locations that are hit too early
      if (config.state.gen < params.startGen) {
        for (unsigned s = 0; s < catalysts.size(); s++) {
          LifeState hitLocations = config.state.Convolve(catalysts[s].reactionMask);
          masks[s].Join(hitLocations);
        }
      }
      if (config.count == params.numCatalysts) {
        bool allRecovered = true;
        for (unsigned i = 0; i < config.count; i++) {
          if (!hasRecovered[i] || missingTime[i] > 0) {
            allRecovered = false;
          }
        }
        if(allRecovered) {
          ReportSolution(config, g, shiftedTargets);
          return;
        }
      }

      // save at least every 15 minutes or so
      if (config.count == 0 || (config.count == 1 && difftime(time(NULL), lastReport) > 900) ){
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
    std::cout << "Usage CatForce.exe <in file>" << std::endl;
    exit(0);
  }

  std::cout << "Input: " << argv[1] << std::endl
            << "Initializing please wait..." << std::endl;

  CatalystSearcher searcher;
  searcher.Init(argv[1]);

  clock_t initialized = clock();
  printf("Total elapsed time: %f seconds\n",
         (double) (initialized - searcher.begin) / CLOCKS_PER_SEC);
  
  std::cout << std::endl
            << "Initialization finished, searching..." << std::endl
            << std::endl;

  searcher.Search();
  // Print report one final time (update files with the final results).
  searcher.Report();
  printf("\n\nFINISH\n");
  clock_t end = clock();
  printf("Total elapsed time: %f seconds\n",
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