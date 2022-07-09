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

const int MAIN_STEP = 1;
__attribute__((flatten)) void MainRun(LifeState &state) {
  state.Step();
}

// const int MAIN_STEP = 2;
// __attribute__((flatten)) void MainChunk(LifeState &state) {
//   state.Step();
// }
// void MainRun(LifeState &state) {
//   MainChunk(state);
//   MainChunk(state);
// }

// const int MAIN_STEP = 4;
// __attribute__((flatten)) void MainChunk(LifeState &state) {
//   state.Step();
//   state.Step();
// }
// void MainRun(LifeState &state) {
//   MainChunk(state);
//   MainChunk(state);
// }

// const int MAIN_STEP = 8;
// __attribute__((flatten)) void MainChunk() {
//   state.Step();
//   state.Step();
// }

// void MainRun() {
//   MainChunk();
//   MainChunk();
//   MainChunk();
//   MainChunk();
// }

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
  int numCatalysts;
  int numTransparent;
  int stableInterval;
  std::string pat;
  int xPat;
  int yPat;
  int startGen;
  int lastGen;
  std::string outputFile;
  std::string fullReportFile;
  int searchArea[4]{};
  int maxW;
  int maxH;
  std::vector<SymmetryTransform> symmetryChain;
  std::vector<std::string> targetFilter;
  std::vector<int> filterdx;
  std::vector<int> filterdy;
  std::vector<int> filterGen;
  std::vector<std::pair<int, int>> filterGenRange;

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
    lastGen = -1;
    outputFile = "results.rle";
    maxW = -1;
    maxH = -1;
    symmetryChain = {};
    maxCatSize = -1;
    fullReportFile = "";
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
  bool mustInclude;
  int period;

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
    period = 1;

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
      } else if (elems[argi] == "mustinclude") {
        mustInclude = true;
        argi += 1;
      } else if (elems[argi] == "period") {
        period = atoi(elems[argi + 1].c_str());
        argi += 2;
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

void CharToTransVec(char ch, std::vector<SymmetryTransform> &trans) {
  trans.push_back(Identity);

  if (ch == '.')
    return;

  if (ch == '|') {
    trans.push_back(ReflectY);
    return;
  }

  if (ch == '-') {
    trans.push_back(ReflectX);
    return;
  }

  if (ch == '+') {
    trans.push_back(ReflectX);
    trans.push_back(ReflectY);
    trans.push_back(Rotate180OddBoth);
    return;
  }

  if (ch == '/' || ch == '\\') {
    trans.push_back(ReflectYeqX);
    return;
  }
  // For 180 degree symetrical
  if (ch == 'x') {
    trans.push_back(Rotate90);
    trans.push_back(ReflectX);
    trans.push_back(ReflectYeqX);
    return;
  }

  if (ch == '*') {
    trans.push_back(ReflectX);
    trans.push_back(ReflectY);
    trans.push_back(Rotate90);
    trans.push_back(Rotate180OddBoth);
    trans.push_back(Rotate270);
    trans.push_back(ReflectYeqX);
    trans.push_back(ReflectYeqNegXP1);
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

  std::string line;

  bool badSymmetry = false;

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
      }

      if (elems[0] == maxWH) {
        params.maxW = atoi(elems[1].c_str());
        params.maxH = atoi(elems[2].c_str());
      }

      if (elems[0] == maxCatSize) {
        params.maxCatSize = atoi(elems[1].c_str());
      }

     std::string symmetryString = "";
      if (elems[0] == symmetry) {
        // reverse-compatibility reasons.
        if (elems[1] == "horizontal") {
          symmetryString = "D2|odd";
        } else if (elems[1] == "horizontaleven") {
          symmetryString = "D2|even";
        } else if (elems[1] == "diagonal") {
          symmetryString = "D2/"; // I think this was the way that it worked before?
        } else if (elems[1] == "rotate180") {
          symmetryString = "C2";
        } else if (elems[1] == "rotate180evenx") {
          symmetryString = "C2horizontaleven";
        } else if (elems[1] == "rotate180evenboth") {
          symmetryString = "C2evenboth";
        } else {
          symmetryString = elems[1];
        }

        std::string start = symmetryString.substr(0,2);
        std::string rest = symmetryString.substr(2);
        if (start == "D2"){
          if (rest == "-" or rest == "vertical" or rest == "verticalodd" or rest == "-odd"){
            params.symmetryChain = {ReflectX};
          } else if (rest == "-even" or rest == "verticaleven"){
            params.symmetryChain = {ReflectXEven};
          } else if (rest == "|" or rest == "horizontal" or rest == "horizontalodd" or rest == "|odd"){
            params.symmetryChain = {ReflectY};
          } else if (rest == "|even" or rest == "horizontaleven"){
            params.symmetryChain = {ReflectYEven};
          } else if ( rest == "/" or rest == "/odd" or rest == "negdiagodd") {
            params.symmetryChain = {ReflectYeqNegX};
          } else if ( rest == "\\" or rest == "\\odd" or rest == "diagodd") {
            params.symmetryChain = {ReflectYeqX};
          } else {
            badSymmetry = true;
          }
        } else if (start == "C2") {
          if (rest == "odd" or rest == "oddboth" or rest == "bothodd" or rest == ""){
            params.symmetryChain = { Rotate180OddBoth};
          } else if (rest == "even" or rest == "botheven" or rest == "evenboth"){
            params.symmetryChain = { Rotate180EvenBoth};
          } else if (rest == "horizontaleven" or rest == "|even"){
            params.symmetryChain = { Rotate180EvenX};
          } else if (rest == "verticaleven" or rest == "-even"){
            params.symmetryChain = {Rotate180EvenY};
          } else {
            badSymmetry = true;
          }
        } else if (start == "C4"){
          if (rest == "" or rest == "odd" or rest == "oddboth" or rest == "bothodd"){
            params.symmetryChain = {Rotate90, Rotate90, Rotate90}; //{Rotate90, Rotate180OddBoth, Rotate270};
          } else if (rest == "even" or rest =="evenboth" or rest == "botheven") {
            params.symmetryChain = { Rotate90Even, Rotate90Even, Rotate90Even};//{ Rotate90Even, Rotate180EvenBoth, Rotate270Even};
          } else {
            badSymmetry = true;
          }
        } else if (start == "D4"){
          std::string evenOddInfo = rest.substr(1);
          if (rest[0] == '+'){
            if(evenOddInfo == "" or evenOddInfo == "odd" or evenOddInfo == "oddboth" or evenOddInfo == "bothodd"){
              params.symmetryChain = { ReflectX, ReflectY, ReflectX};//{ ReflectX, ReflectY, Rotate180OddBoth};
            } else if (evenOddInfo == "even" or evenOddInfo =="evenboth" or evenOddInfo == "botheven"){
              params.symmetryChain = { ReflectXEven, ReflectYEven, ReflectXEven}; //{ ReflectXEven, ReflectYEven, Rotate180EvenBoth};
            } else if ( evenOddInfo == "verticaleven" or evenOddInfo == "-even") {
              params.symmetryChain = { ReflectXEven, ReflectY, ReflectXEven};//{ ReflectXEven, ReflectY, Rotate180EvenX};
            } else if ( evenOddInfo == "horizontaleven" or evenOddInfo == "|even") {
              params.symmetryChain = { ReflectX, ReflectYEven, ReflectX};//{ ReflectX, ReflectYEven, Rotate180EvenY};
            } else {
              badSymmetry = true;
            }
          } else if (rest[0] == 'x') {
            if (evenOddInfo == "odd" or evenOddInfo == "oddboth" or evenOddInfo == ""){
              params.symmetryChain = {ReflectYeqX, ReflectYeqNegXP1, ReflectYeqX};//{ ReflectYeqX, ReflectYeqNegXP1, Rotate180OddBoth};
            } else if (evenOddInfo == "even" or evenOddInfo == "evenboth"){
              params.symmetryChain = { ReflectYeqX, ReflectYeqNegX,ReflectYeqX};//{ ReflectYeqX, ReflectYeqNegX, Rotate180EvenBoth};
            } else {
              badSymmetry = true;
            }
          } else {
              badSymmetry = true;
          }
        } else if (start == "D8") {
          if (rest == "odd" or rest == "oddboth" or rest == ""){
            params.symmetryChain = {ReflectY, ReflectYeqNegXP1, ReflectX, ReflectYeqX, ReflectY, ReflectYeqNegXP1, ReflectX};
            // reflections are faster than 90 degree rotations, so we reflect around in a circle.
            //{ ReflectX, ReflectY, Rotate90, Rotate270, Rotate180OddBoth, ReflectYeqX, ReflectYeqNegXP1};
          } else if (rest == "even" or rest == "evenboth"){
            params.symmetryChain = {ReflectYEven, ReflectYeqNegX, ReflectXEven, ReflectYeqX, ReflectYEven, ReflectYeqNegX, ReflectXEven};
            //{ ReflectXEven, ReflectYEven, Rotate90Even, Rotate270Even, Rotate180EvenBoth, ReflectYeqX, ReflectYeqNegX};
          } else {
            badSymmetry = true;
          }
        } else {
          badSymmetry = true;
        }
      }

    } catch (const std::exception &ex) {
    }
  }
  if(params.lastGen == -1)
    params.lastGen = params.maxGen - 1;

  if (params.pat.length() == 0) {
    std::cout << "Did not read any pattern!" << std::endl;
    exit(1);
  }
  if (catalysts.empty()) {
    std::cout << "Did not read any catalysts!" << std::endl;
    exit(1);
  }
  if (badSymmetry) {
    std::cout << "Couldn\'t parse symmetry option" << std::endl;
    exit(1);
  }
}

class CatalystData {
public:
  LifeState state;
  std::vector<LifeState> phases;
  // LifeTarget target;
  LifeState reactionMask;
  std::vector<LifeState> phaseReactionMask;
  int maxDisappear;
  std::vector<LifeTarget> forbidden;
  LifeState required;
  LifeState locus;
  bool transparent;
  bool mustInclude;
  int period;

  static std::vector<CatalystData> FromInput(CatalystInput &input);
};

std::vector<CatalystData> CatalystData::FromInput(CatalystInput &input) {
  std::vector<SymmetryTransform> trans;
  CharToTransVec(input.symmType, trans);

  const char *rle = input.rle.c_str();

  std::vector<CatalystData> results;

  for (auto &tran : trans) {
    LifeState pat = LifeState::Parse(rle, input.centerX, input.centerY, tran);

    for (int i = 0; i < input.period; i++) {
      CatalystData result;

      result.state = pat;
      result.reactionMask = pat.BigZOI();
      result.reactionMask.Transform(Rotate180OddBoth);

      // result.target = LifeTarget(pat);
      result.maxDisappear = input.maxDisappear;

      for (int k = 0; k < input.forbiddenRLE.size(); k++) {
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

      LifeState tmp = pat;
      for (int j = 0; j < input.period; j++) {
        result.phases.push_back(tmp);

        LifeState phasemask = tmp.BigZOI();
        phasemask.Transform(Rotate180OddBoth);
        result.phaseReactionMask.push_back(phasemask);

        tmp.Step();
      }

      result.transparent = input.transparent;
      result.period = input.period;
      result.mustInclude = input.mustInclude;

      results.push_back(result);

      pat.Step();
    }
  }
  return results;
}

struct Configuration {
  int count;
  int transparentCount;
  int mustIncludeCount;
  std::array<int, MAX_CATALYSTS> curx;
  std::array<int, MAX_CATALYSTS> cury;
  std::array<int, MAX_CATALYSTS> curs;
  // int minIter;
  LifeState state;
  LifeState catalystsState;
  LifeState startingCatalysts;
  // std::vector<LifeTarget> shiftedTargets;
};

// Fix a, what positions of b causes a collision?
LifeState CollisionMask(const LifeState &a, const LifeState &b, int period) {
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

      for(int i = 0; i < period; i++)
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
    LifeState mask = CollisionMask(a.state, b.state, a.period * b.period);
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

    // for (int param : params)
    //   std::cout << param << ",";

    std::cout << std::endl;
  }
};

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

    int howmany = results.size();

    if (maxCatSize != -1)
      howmany = std::min(howmany, maxCatSize);

    int width = Dist * howmany;
    int height = Dist;

    std::vector<std::vector<int>> vec(width, std::vector<int>(height));

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

  void Add(LifeState &init, const LifeState &afterCatalyst, const LifeState &catalysts,
           const Configuration &conf, int firstGenSurvive,
           int genSurvive) {
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
  std::vector<CatalystData> catalysts;
  std::vector<LifeTarget> targetFilter;
  std::vector<std::vector<LifeState>> catalystCollisionMasks;

  clock_t current{};
  long long idx{};
  int found{};
  int fullfound{};
  long long total{};
  unsigned short int counter{};

  CategoryContainer *categoryContainer{};
  CategoryContainer *fullCategoryContainer{};

  bool hasFilter{};
  bool reportAll{};
  bool hasFilterDontReportAll{};

  int filterMaxGen{};

  void Init(const char *inputFile, int nthreads) {
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

    for (int i = 0; i < params.targetFilter.size(); i++)
      targetFilter.push_back(LifeTarget::Parse(params.targetFilter[i].c_str(),
                                               params.filterdx[i], params.filterdy[i]));

    catalystCollisionMasks = std::vector<std::vector<LifeState>>(
        catalysts.size(), std::vector<LifeState>(catalysts.size()));

    for (int s = 0; s < catalysts.size(); s++) {
      // LifeState nonLocus = catalysts[s];
      // nonLocus.Copy(catalystLocus[s], ANDNOT);

      // catalystAvoidMasks[s] = nonLocus.BigZOI();
      // catalystAvoidMasks[s].Transform(Rotate180OddBoth);

      // catalystLocusReactionMasks[s] = catalystLocus[s].BigZOI();
      // catalystLocusReactionMasks[s].Transform(Rotate180OddBoth);
      // catalystLocusReactionMasks[s].Copy(catalystAvoidMasks[s], ANDNOT);

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

    int fact = 1;

    for (int i = 0; i < catalysts.size(); i++) {
      total *= params.searchArea[2];
      total *= params.searchArea[3];
      total *= catalysts.size();
      fact *= (i + 1);
    }

    total /= fact;

    std::cout << "Approximated Total: " << total << std::endl;
    total = total / 1000000;

    if (total == 0)
      total++;

    hasFilter = !params.targetFilter.empty();
    reportAll = params.fullReportFile.length() != 0;

    filterMaxGen = FilterMaxGen();
  }

  int FilterMaxGen() {
    int maxGen = -1;

    for (int j = 0; j < targetFilter.size(); j++) {
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

  void IncreaseIndexAndReport(bool saveFile = true) {
    counter++;

    if (counter == 0) {
      idx += 65536;

      if (idx % (1048576) == 0) {
        if ((double) (clock() - current) / CLOCKS_PER_SEC > 10) {
          current = clock();
          Report(saveFile);
        }
      }
    }
  }

  bool HasForbidden(Configuration &conf, int curIter) {
    LifeState workspace;
    workspace.Join(conf.startingCatalysts);
    workspace.JoinWSymChain(pat, params.symmetryChain);

    for (int i = 0; i <= curIter + 1; i++) {
      for (int j = 0; j < catalysts.size(); j++) {
        for (int k = 0; k < catalysts[conf.curs[j]].forbidden.size(); k++) {
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
    for (int j = 0; j < targetFilter.size(); j++) {
      if (workspace.gen == params.filterGen[j] &&
          workspace.Contains(targetFilter[j]) == false) {
        return true;
      }
    }

    return false;
  }

  bool ValidateFilters(Configuration &conf) {
    LifeState workspace;
    workspace.JoinWSymChain(pat, params.symmetryChain);
    workspace.Join(conf.startingCatalysts);

    std::vector<bool> rangeValid(params.filterGen.size(), false);

    for (int k = 0; k < params.filterGen.size(); k++)
      if (params.filterGen[k] >= 0)
        rangeValid[k] = true;

    for (int j = 0; j <= filterMaxGen; j++) {
      for (int k = 0; k < params.filterGen.size(); k++) {
        if (workspace.gen == params.filterGen[k] &&
            workspace.Contains(targetFilter[k]) == false)
          return false;

        if (params.filterGen[k] == -1 &&
            params.filterGenRange[k].first <= workspace.gen &&
            params.filterGenRange[k].second >= workspace.gen &&
            workspace.Contains(targetFilter[k]) == true)
          rangeValid[k] = true;
      }

      workspace.Step();
    }

    for (int k = 0; k < params.filterGen.size(); k++)
      if (!rangeValid[k])
        return false;

    return true;
  }

  void ReportSolution(Configuration &conf, int successtime){
    LifeState init;
    LifeState afterCatalyst;
    LifeState catalysts;

    LifeState workspace;
    // if reportAll - ignore filters and update fullReport
    if (reportAll) {
      workspace.Join(conf.startingCatalysts);
      workspace.JoinWSymChain(pat, params.symmetryChain);
      init.Copy(workspace);

      workspace.Step(successtime - params.stableInterval + 2);
      afterCatalyst.Copy(workspace);

      fullfound++;

      fullCategoryContainer->Add(init, afterCatalyst, conf.startingCatalysts, conf,
                                 successtime - params.stableInterval + 2, 0);
    }

    // If has fitlter validate them;
    if (hasFilter) {
      if (!ValidateFilters(conf))
        return;
    }

    // if (HasForbidden(conf, successtime + 3))
    //   return;

    // If all filters validated update results
    workspace.Clear();
    workspace.JoinWSymChain(pat, params.symmetryChain);
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
    config.count = 0;
    config.transparentCount = 0;
    config.mustIncludeCount = 0;
    config.state.JoinWSymChain(pat, params.symmetryChain);
    LifeState history = config.state;

    LifeState bounds =
        LifeState::SolidRect(params.searchArea[0], params.searchArea[1],
                             params.searchArea[2], params.searchArea[3]);

    std::vector<LifeState> masks(catalysts.size());
    for (int s = 0; s < catalysts.size(); s++) {
      masks[s] = config.state.Convolve(catalysts[s].reactionMask);
      masks[s].Copy(bounds, ORNOT);
    }

    std::vector<LifeTarget> shiftedTargets(params.numCatalysts);

    RecursiveSearch(config, history, LifeState(), masks, shiftedTargets,
                    std::array<int, MAX_CATALYSTS>(), std::array<int, MAX_CATALYSTS>(),
                    std::array<bool, MAX_CATALYSTS>(), std::array<bool, MAX_CATALYSTS>());
  }

  void
  RecursiveSearch(Configuration config, LifeState history, const LifeState required,
                  std::vector<LifeState> masks,
                  std::vector<LifeTarget> &shiftedTargets, // This can be shared

                  std::array<int, MAX_CATALYSTS> missingTime,
                  std::array<int, MAX_CATALYSTS> recoveredTime,
                  std::array<bool, MAX_CATALYSTS> hasReacted,
                  std::array<bool, MAX_CATALYSTS> hasRecovered) {

    for (int g = config.state.gen; g < params.maxGen; g++) {
      if (config.count == 0 && g > params.lastGen)
        return;

      if (config.count == 0) {
        std::cout << "Collision at gen " << g << std::endl;
      }

      if (!config.state.Contains(required))
        return;

      for (int i = 0; i < config.count; i++) {
        // if (hasRecovered[i]) {
        //   continue;
        // }
        if (config.state.gen % catalysts[config.curs[i]].period != 0)
          continue;

        if (config.state.Contains(shiftedTargets[i])) {
          missingTime[i] = 0;
          recoveredTime[i] += catalysts[config.curs[i]].period;
        } else {
          hasReacted[i] = true;
          missingTime[i] += catalysts[config.curs[i]].period;
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
        newcells.Copy(config.catalystsState, ANDNOT);

        if (!newcells.IsEmpty()) {
          // for (int s = 0; s < catalysts.size(); s++) {
          //   LifeState hitLocations = newcells.Convolve(catalystAvoidMasks[s]);
          //   masks[s].Join(hitLocations);
          // }

          for (int s = 0; s < catalysts.size(); s++) {
            if (config.transparentCount == params.numTransparent && catalysts[s].transparent)
              continue;
            if (config.count == params.numCatalysts - 1 && config.mustIncludeCount == 0 && !catalysts[s].mustInclude)
              continue;

            LifeState newPlacements =
                catalysts[s]
                    .phaseReactionMask[g % catalysts[s].period]
                    .Convolve(newcells);
            newPlacements.Copy(masks[s], ANDNOT);

            while (!newPlacements.IsEmpty()) {
              // Do the placement
              auto newPlacement = newPlacements.FirstOn();

              if (config.count == 0) {
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
              LifeState symCatalyst;
              symCatalyst.JoinWSymChain(shiftedCatalyst, params.symmetryChain);
              newConfig.startingCatalysts.Join(symCatalyst);
              for (int k = 0; k < g % catalysts[s].period; k++)
                symCatalyst.Step();
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

              shiftedTargets[config.count] = LifeTarget(shiftedCatalyst);

              LifeState bounds;
              if (params.maxW != -1) {
                bounds = LifeState::SolidRect(newPlacement.first - params.maxW,
                                              newPlacement.second - params.maxH,
                                              2 * params.maxW - 1,
                                              2 * params.maxH - 1);
              }

              // If we just placed the last catalyst, don't bother
              if (newConfig.count != params.numCatalysts) {
                for (int t = 0; t < catalysts.size(); t++) {
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
        for (int s = 0; s < catalysts.size(); s++) {
          LifeState hitLocations = config.state.Convolve(catalysts[s].phaseReactionMask[g % catalysts[s].period]);
          masks[s].Join(hitLocations);
        }
      }
      if (config.count == params.numCatalysts) {
        bool allRecovered = true;
        for (int i = 0; i < config.count; i++) {
          if (!hasRecovered[i] || missingTime[i] > 0) {
            allRecovered = false;
          }
        }
        if(allRecovered) {
          ReportSolution(config, g);
          return;
        }
      }

      if (config.count == 0)
        Report();

      history.Copy(config.state, OR);
      config.state.Step();
      config.catalystsState.Step();
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
}
