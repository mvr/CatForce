// CatForce - Catalyst search utility based on LifeAPI using brute force.
// Written by Michael Simkin 2015
#include <omp.h>
#include "LifeAPI.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <ctime>
#include <vector>


const int none[] = {1, 0, 0, 1};
const int flipX[] = {-1, 0, 0, 1};
const int flipY[] = {1, 0, 0, -1};
const int flipXY[] = {-1, 0, 0, -1};

const int rot90clock[] = {0, 1, -1, 0};
const int rot90anti[] = {0, -1, 1, 0};
const int symmXY[] = {0, 1, 1, 0};
const int symmYX[] = {0, -1, -1, 0};

const int MAIN_STEP = 1;

__attribute__((flatten)) void MainRun(LifeState &state) {
  IterateState(&state);
}

// const int MAIN_STEP = 2;
// __attribute__((flatten)) void MainRun() {
//   IterateState(GlobalState);
//   IterateState(GlobalState);
// }

// const int MAIN_STEP = 4;
// __attribute__((flatten)) void MainChunk() {
//   IterateState(GlobalState);
//   IterateState(GlobalState);
// }
// void MainRun() {
//   MainChunk();
//   MainChunk();
// }

// const int MAIN_STEP = 8;
// __attribute__((flatten)) void MainChunk() {
//   IterateState(GlobalState);
//   IterateState(GlobalState);
// }

// void MainRun() {
//   MainChunk();
//   MainChunk();
//   MainChunk();
//   MainChunk();
// }

enum Symmetry {
  NONE,
  HORIZONTAL,
  HORIZONTALEVEN,
  DIAGONAL,
  ROTATE180,
  ROTATE180EVENX,
  ROTATE180EVENBOTH
};

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
  Symmetry symmetricSearch;
  std::vector<std::string> targetFilter;
  std::vector<int> filterdx;
  std::vector<int> filterdy;
  std::vector<int> filterGen;
  std::vector<std::pair<int, int>> filterGenRange;

  int maxCatSize;
  bool combineResults;
  std::vector<int> combineSurvive;

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
    lastGen = 100000;
    outputFile = "results.rle";
    maxW = -1;
    maxH = -1;
    symmetricSearch = NONE;
    maxCatSize = -1;
    fullReportFile = "";
    combineResults = false;
  }
};

class CatalystInput {
public:
  std::string rle;
  int maxDesapear;
  int centerX;
  int centerY;
  char symmType;
  std::vector<std::string> forbiddenRLE;
  std::vector<std::pair<int, int>> forbiddenXY;

  explicit CatalystInput(std::string &line) {
    std::vector<std::string> elems;
    split(line, ' ', elems);

    if (elems.size() < 6) {
      std::cout << "The line " << line << "is invalid" << std::endl;
      std::cout << "Format: cat <rle> <absense interval> <centerX> <centerY> "
                   "<symm Type | + / x *>"
                << std::endl;
      getchar();
      exit(0);
    }

    rle = elems[1];
    maxDesapear = atoi(elems[2].c_str());
    centerX = atoi(elems[3].c_str());
    centerY = atoi(elems[4].c_str());
    symmType = elems[5].at(0);

    int argi = 6;

    while (argi + 3 < elems.size()) {
      if (elems[argi] == "forbidden") {
        forbiddenRLE.push_back(elems[argi + 1]);
        forbiddenXY.emplace_back(
            atoi(elems[argi + 2].c_str()), atoi(elems[argi + 3].c_str()));

        argi += 4;
      }
    }
  }

  void Print() const {
    std::cout << rle << " " << maxDesapear << " " << centerX << " " << centerY
              << " " << symmType << std::endl;
  }
};

void CharToTransVec(char ch, std::vector<const int *> &trans) {
  trans.push_back(none);

  if (ch == '.')
    return;

  if (ch == '|') {
    trans.push_back(flipY);
    return;
  }

  if (ch == '-') {
    trans.push_back(flipX);
    return;
  }

  if (ch == '+') {
    trans.push_back(flipX);
    trans.push_back(flipY);
    trans.push_back(flipXY);
    return;
  }

  if (ch == '/' || ch == '\\') {
    trans.push_back(symmXY);
    return;
  }
  // For 180 degree symetrical
  if (ch == 'x') {
    trans.push_back(rot90clock);
    trans.push_back(flipX);
    trans.push_back(symmXY);
    return;
  }

  if (ch == '*') {
    trans.push_back(flipX);
    trans.push_back(flipY);
    trans.push_back(flipXY);
    trans.push_back(symmYX);
    trans.push_back(symmXY);
    trans.push_back(rot90anti);
    trans.push_back(rot90clock);
    return;
  }
}

void ReadParams(const std::string& fname, std::vector<CatalystInput> &catalysts,
                SearchParams &params) {
  std::ifstream infile;
  infile.open(fname.c_str(), std::ifstream::in);
  if (!infile.good()) {
    std::cout << "Could not open file!" << std::endl;
    exit(0);
  }

  std::string Cat = "cat";
  std::string maxGen = "max-gen";
  std::string startGen = "start-gen";
  std::string lastGen = "last-gen";

  std::string numCat = "num-catalyst";
  std::string stable = "stable-interval";
  std::string area = "search-area";
  std::string pat = "pat";
  std::string outputFile = "output";
  std::string filter = "filter";
  std::string maxWH = "fit-in-width-height";
  std::string maxCatSize = "max-category-size";
  std::string fullReport = "full-report";
  std::string combine = "combine-results";

  std::string symmetry = "symmetry";

  std::string line;

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

      if (elems[0] == combine && elems[1] == "yes") {
        params.combineResults = true;

        for (int i = 2; i < elems.size(); i++)
          params.combineSurvive.push_back(atoi(elems[i].c_str()));
      }

      if (elems[0] == symmetry) {
        if (elems[1] == "horizontal") {
          params.symmetricSearch = HORIZONTAL;
        } else if (elems[1] == "horizontaleven") {
          params.symmetricSearch = HORIZONTALEVEN;
        } else if (elems[1] == "diagonal") {
          params.symmetricSearch = DIAGONAL;
        } else if (elems[1] == "rotate180") {
          params.symmetricSearch = ROTATE180;
        } else if (elems[1] == "rotate180evenx") {
          params.symmetricSearch = ROTATE180EVENX;
        } else if (elems[1] == "rotate180evenboth") {
          params.symmetricSearch = ROTATE180EVENBOTH;
        }
      }

    } catch (const std::exception &ex) {
    }
  }
  if (params.pat.length() == 0) {
    std::cout << "Did not read any pattern!" << std::endl;
    exit(0);
  }
  if (catalysts.empty()) {
    std::cout << "Did not read any catalysts!" << std::endl;
    exit(0);
  }
}

void ApplySym(LifeState *state, Symmetry sym) {
  if (sym == HORIZONTAL) {
    FlipX(state);
  } else if (sym == HORIZONTALEVEN) {
    Reverse(state->state, 0, N - 1);
  } else if (sym == ROTATE180) { // before: DIAGONAL
    FlipX(state);
    FlipY(state);
  } else if (sym == ROTATE180EVENX) { // before: DIAGONALEVENX
    Reverse(state->state, 0, N - 1);
    FlipY(state);
  } else if (sym == ROTATE180EVENBOTH) { // before: DIAGONALEVENBOTH
    Reverse(state->state, 0, N - 1);
    BitReverse(state);
  } else if (sym == DIAGONAL) {
    Transpose(state);
  }
}

void PutStateWSym(LifeState *main, LifeState *state, int x, int y, Symmetry sym) {
  Join(main, state, x, y);
  if (sym == NONE)
    return;
  LifeState transformed;
  ClearData(&transformed);

  transformed.min = 0;
  transformed.max = N - 1;
  transformed.gen = 0;

  Join(&transformed, state, x, y);
  ApplySym(&transformed, sym);
  Join(main, &transformed);
}

void PutStateWSym(LifeState *main, const LifeState *state, Symmetry sym) {
  Join(main, state);
  if (sym == NONE)
    return;

  LifeState transformed;
  Copy(&transformed, state);
  ApplySym(&transformed, sym);
  Join(main, &transformed);
}

void GenerateStates(const std::vector<CatalystInput> &catalysts,
                    std::vector<LifeState *> &states,
                    std::vector<std::vector<LifeTarget *>> &forbidden,
                    std::vector<int> &maxSurvive) {
  for (const auto & catalyst : catalysts) {
    std::vector<const int *> trans;
    CharToTransVec(catalyst.symmType, trans);

    const char *rle = catalyst.rle.c_str();
    int dx = catalyst.centerX;
    int dy = catalyst.centerY;
    int maxDesapear = catalyst.maxDesapear;

    for (auto & tran : trans) {
      int dxx = tran[0];
      int dxy = tran[1];
      int dyx = tran[2];
      int dyy = tran[3];

      states.push_back(NewState(rle, dx, dy, dxx, dxy, dyx, dyy));
      maxSurvive.push_back(maxDesapear);

      std::vector<LifeTarget *> forbidTarg;

      for (int k = 0; k < catalyst.forbiddenRLE.size(); k++) {
        forbidTarg.push_back(NewTarget(
            NewState(catalyst.forbiddenRLE[k].c_str(),
                     catalyst.forbiddenXY[k].first,
                     catalyst.forbiddenXY[k].second, dxx, dxy, dyx, dyy)));
      }

      forbidden.push_back(forbidTarg);
    }
  }
}

void InitCatalysts(const std::string& fname, std::vector<LifeState *> &states,
                   std::vector<std::vector<LifeTarget *>> &forbidden,
                   std::vector<int> &maxSurvive, SearchParams &params) {
  std::vector<CatalystInput> catalysts;
  ReadParams(fname, catalysts, params);
  GenerateStates(catalysts, states, forbidden, maxSurvive);
}

void XYStartGenPerState(const std::vector<LifeTarget *> &targets,
                        const LifeState *pat, const SearchParams &params,
                        const std::vector<LifeState *> &states,
                        std::vector<std::vector<std::vector<int>>> &statexyGen,
                        int const nthreads) {
  if (params.numCatalysts == 1) {
    for (long i = 0; i < states.size(); i++) {
      std::vector<std::vector<int>> xyVec;
      xyVec.reserve(64);

      for (int x = 0; x < 64; x++) {
        std::vector<int> xVec(64, params.startGen);
        xyVec.push_back(xVec);
      }
      statexyGen.push_back(xyVec);
    }
    return;
  }


  const int chunksize = states.size() / ((unsigned long)nthreads);
  std::vector<std::pair<long, long>> chunkbounds;
  long lowerbound = 0;
  for (int chunki = 0; chunki < nthreads-1; chunki++){
      chunkbounds.emplace_back(lowerbound, lowerbound + chunksize);
      lowerbound += chunksize;
  }
  chunkbounds.emplace_back(lowerbound, states.size());
  statexyGen.reserve(states.size());
  #pragma omp parallel for ordered schedule(static,1) default(none) shared(states, params, targets, statexyGen, pat, chunkbounds, std::cout)
  for (auto & bounds : chunkbounds){
    std::vector<std::vector<std::vector<int>>> perthread_statexyGen;
    perthread_statexyGen.reserve(bounds.second-bounds.first);

    #pragma omp critical
    for (long i = bounds.first; i < bounds.second; i++) {
      std::vector<std::vector<int>> xyVec;
      xyVec.reserve(64);

      for (int x = 0; x < 64; x++) {
        std::vector<int> xVec;
        xVec.reserve(64);

        for (int y = 0; y < 64; y++) {
          LifeState state;
          ClearData(&state);
          PutStateWSym(&state, states[i], x, y, params.symmetricSearch);
          PutStateWSym(&state, pat, params.symmetricSearch);
          int j;

          for (j = 0; j < params.maxGen + 5; j++) {
            if (Contains(&state, targets[i], x, y) == NO) {
              break;
            }
            Run(&state, 1);
          }

          if (j == params.maxGen + 4)
            j = -1;

          xVec.push_back(j - 1);
        }
        xyVec.push_back(xVec);
      }

      perthread_statexyGen.push_back(xyVec);
    }

    #pragma omp ordered
    {
      statexyGen.insert(
          statexyGen.end(),
          std::make_move_iterator(perthread_statexyGen.begin()),
          std::make_move_iterator(perthread_statexyGen.end())
      );
    };
  }
}

void PreIteratePat(LifeState *pat, std::vector<LifeState *> &preIterated,
                   const SearchParams &params) {
  LifeState workspace;
  ClearData(&workspace);
  PutStateWSym(&workspace, pat, params.symmetricSearch);

  for (int i = 0; i < params.maxGen + 5; i++) {
    LifeState *t = NewState();
    Copy(t, &workspace);
    preIterated.push_back(t);
    Run(&workspace, 1);
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

class SearchResult {
public:
  // Saved for the report
  LifeState *init;

  // iters state in form of integers
  std::vector<int> params;
  int maxGenSurvive;
  int firstGenSurvive;

  SearchResult(LifeState *initState, const Configuration &conf,
               int firstGenSurviveIn, int genSurvive) {
    init = NewState();
    Copy(init, initState);

    for (int i = 0; i < conf.curs.size(); i++) {
      params.push_back(conf.curs[i]);
      params.push_back(conf.curx[i]);
      params.push_back(conf.cury[i]);
    }

    maxGenSurvive = genSurvive;
    firstGenSurvive = firstGenSurviveIn;
  }

  // int SetIters(Enumerator &enu, const int &startIdx) {
  //   int idx = startIdx;

  //   for (int i = 0; i < params.size(); i += 3) {
  //     enu.curs[idx] = params[i];
  //     enu.curx[idx] = params[i + 1];
  //     enu.cury[idx] = params[i + 2];
  //     idx++;
  //   }

  //   return idx;
  // }

  void Print() {
    std::cout << "start:" << firstGenSurvive;
    std::cout << ", finish:" << maxGenSurvive << ", params: ";

    for (int param : params)
      std::cout << param << ",";

    std::cout << std::endl;
  }

  ~SearchResult() { FreeState(init); }
};

class Category {
private:
  LifeState *tempCat;
  LifeState *tempTest;
  int catDelta;
  int maxgen;
  uint64_t hash;

public:
  LifeState *categoryKey;
  std::vector<SearchResult *> results;

  Category(LifeState *catalystRemoved, SearchResult *firstResult,
           int catDeltaIn, int maxGen) {
    categoryKey = NewState();
    Copy(categoryKey, catalystRemoved);
    results.push_back(firstResult);
    catDelta = catDeltaIn;
    tempCat = NewState();
    tempTest = NewState();
    maxgen = maxGen;

    ClearData(tempCat);
    Copy(tempCat, categoryKey);
    Evolve(tempCat, maxgen - tempCat->gen);
    hash = GetHash(tempCat);
  }

  void Add(SearchResult *result) {
    if(results.size() < 100)
      results.push_back(result);
  }

  bool BelongsTo(LifeState *test, const uint64_t &testHash) {
    if (testHash != hash)
      return false;

    ClearData(tempCat);
    Copy(tempCat, categoryKey);

    ClearData(tempTest);
    Copy(tempTest, test);

    if (tempTest->gen > tempCat->gen)
      Evolve(tempCat, tempTest->gen - tempCat->gen);
    else if (tempTest->gen < tempCat->gen)
      Evolve(tempTest, tempCat->gen - tempTest->gen);

    for (int i = 0; i < catDelta; i++) {
      if (AreEqual(tempTest, tempCat) == YES)
        return true;

      Evolve(tempTest, 1);
      Evolve(tempCat, 1);
    }

    return false;
  }

  static bool CompareSearchResult(SearchResult *a, SearchResult *b) {
    return (a->maxGenSurvive - a->firstGenSurvive) >
           (b->maxGenSurvive - b->firstGenSurvive);
  }

  void Sort() {
    std::sort(results.begin(), results.end(), Category::CompareSearchResult);
  }

  void RemoveTail() {
    if (results.size() <= 1)
      return;

    for (int i = 1; i < results.size(); i++)
      delete results[i];

    results.erase(results.begin() + 1, results.end());
  }

  ~Category() {
    FreeState(tempCat);
    FreeState(tempTest);

    for (auto & result : results)
      delete result;
  }

  void Print() {
    for (auto & result : results)
      result->Print();
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
          vec[Dist * l + i][j] = Get(i, j, results[l]->init->state);

    return GetRLE(vec);
  }
};

class CategoryContainer {
public:
  std::vector<Category *> categories;
  LifeState *tempState;
  int catDelta;
  int maxgen;

  explicit CategoryContainer(int maxGen) {
    catDelta = 14;
    tempState = NewState();
    maxgen = maxGen + catDelta;
  }

  CategoryContainer(int cats, int maxGen) {
    catDelta = cats;
    tempState = NewState();
    maxgen = maxGen + catDelta;
  }

  void Add(LifeState *init, LifeState *afterCatalyst, LifeState *catalysts,
           const Configuration &conf, int firstGenSurvive,
           int genSurvive) {
    LifeState *result = NewState();

    ClearData(result);
    Copy(result, afterCatalyst);
    Copy(result, catalysts, XOR);

    Evolve(result, maxgen - result->gen);
    uint64_t hash = GetHash(result);

    ClearData(result);
    Copy(result, afterCatalyst);
    Copy(result, catalysts, XOR);

    for (auto & category: categories) {
      if (category->BelongsTo(result, hash)) {
        if (category->results[0]->params.size() == 3 * conf.curs.size())
          category->Add(
              new SearchResult(init, conf, firstGenSurvive, genSurvive));
        return;
      }
    }

    LifeState *categoryKey = NewState();
    Copy(categoryKey, afterCatalyst);
    Copy(categoryKey, catalysts, XOR);

    categories.push_back(new Category(
        categoryKey, new SearchResult(init, conf, firstGenSurvive, genSurvive),
        catDelta, maxgen));
  }

  void Sort() {
    for (auto & category: categories)
      category->Sort();
  }

  void Print() {
    for (auto & category: categories)
      category->Print();
  }

  void RemoveTail() {
    for (auto & category: categories)
      category->RemoveTail();
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
  std::vector<LifeState *> states;
  std::vector<int> maxSurvive;
  SearchParams params;
  LifeState *pat{};
  int numIters{};
  Enumerator enu;
  std::vector<LifeTarget *> targetFilter;
  std::vector<LifeTarget *> targets;
  std::vector<std::vector<LifeTarget *>> forbiddenTargets;
  std::vector<std::vector<std::vector<int>>> statexyGen;
  std::vector<LifeState *> preIterated;
  std::vector<int> activated;
  std::vector<int> absentCount;
  clock_t current{};
  long long idx{};
  int found{};
  int fullfound{};
  long long total{};
  unsigned short int counter{};
  CategoryContainer *categoryContainer{};
  CategoryContainer *fullCategoryContainer{};

  // flags and memeber for the search

  bool hasFilter{};
  bool reportAll{};
  bool hasFilterDontReportAll{};

  int filterMaxGen{};
  int iterationMaxGen{};

  int surviveCountForUpdate{};

  LifeState *init{};
  LifeState *afterCatalyst{};
  LifeState *catalysts{};

  void Init(const char *inputFile, int nthreads) {
    begin = clock();
    InitCatalysts(inputFile, states, forbiddenTargets, maxSurvive, params);
    pat = NewState(params.pat.c_str(), params.xPat, params.yPat);
    numIters = params.numCatalysts;
    categoryContainer = new CategoryContainer(params.maxGen);
    fullCategoryContainer = new CategoryContainer(params.maxGen);

    enu.count = params.numCatalysts;
    enu.x = params.searchArea[0];
    enu.y = params.searchArea[1];
    enu.w = params.searchArea[2];
    enu.h = params.searchArea[3];
    enu.maxW = params.maxW;
    enu.maxH = params.maxH;
    enu.s = states.size();
    enu.states = std::vector<LifeState *>();
    for (int i = 0; i < enu.s; i++) {
      LifeState *state = NewState();
      Copy(state, states[i]);
      enu.states.push_back(state);
    }
    Reset(enu);

    for (int i = 0; i < params.targetFilter.size(); i++)
      targetFilter.push_back(NewTarget(params.targetFilter[i].c_str(),
                                       params.filterdx[i], params.filterdy[i]));

    for (auto & state : states)
      targets.push_back(NewTarget(state));

    XYStartGenPerState(targets, pat, params, states, statexyGen, nthreads);
    enu.activations = statexyGen;
    enu.cumulActivation = std::vector<int>(enu.count, enu.activations[0][(enu.x + 64) % 64][(enu.y + 64) % 64]);
    for (int i = 0; i < enu.s; i++) {
      for (int x = 0; x < 64; x++) {
        for (int y = 0; y < 64; y++) {
          int minIter = statexyGen[i][x][y];
          enu.quickEnough[i][x][y] = minIter < params.lastGen;
          enu.slowEnough[i][x][y] = params.startGen <= minIter;
        }
      }
    }

    PreIteratePat(pat, preIterated, params);
    AddIterators(numIters);

    current = clock();
    idx = 0;
    found = 0;
    total = 1;
    counter = 0;

    fullfound = 0;

    int fact = 1;

    for (int i = 0; i < numIters; i++) {
      total *= params.searchArea[2];
      total *= params.searchArea[3];
      total *= states.size();
      fact *= (i + 1);
    }

    total /= fact;

    std::cout << "Approximated Total: " << total << std::endl;
    total = total / 1000000;

    if (total == 0)
      total++;

    hasFilter = !params.targetFilter.empty();
    reportAll = params.fullReportFile.length() != 0;
    hasFilterDontReportAll = hasFilter && !reportAll;

    filterMaxGen = FilterMaxGen();
    iterationMaxGen = params.maxGen;

    init = NewState();
    afterCatalyst = NewState();
    catalysts = NewState();

    surviveCountForUpdate = params.stableInterval;

    if (params.combineResults) {
      if (params.combineSurvive.empty())
        params.combineSurvive.push_back(1);

      surviveCountForUpdate = params.combineSurvive[0];
      hasFilter = false;
      hasFilterDontReportAll = false;
    }
  }

  void AddIterators(int num) {
    for (int i = 0; i < num; i++) {
      activated.push_back(0);
      absentCount.push_back(0);
    }

    numIters = num;
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
    float percent = (idx / 10000) / (total * 1.0);
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
      std::cout << "Saving " << params.outputFile << std::endl;

      std::ofstream catResultsFile(params.outputFile.c_str());
      catResultsFile << "x = 0, y = 0, rule = B3/S23\n";
      catResultsFile << categoryContainer->CategoriesRLE(params.maxCatSize);
      catResultsFile.close();

      if (params.fullReportFile.length() != 0) {
        std::cout << "Saving " << params.fullReportFile << std::endl;

        std::ofstream fullCatResultsFile(params.fullReportFile.c_str());
        fullCatResultsFile << "x = 0, y = 0, rule = B3/S23\n";
        fullCatResultsFile << fullCategoryContainer->CategoriesRLE(params.maxCatSize);
        fullCatResultsFile.close();
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

  void InitActivationCounters() {
    for (int i = 0; i < numIters; i++) {
      activated[i] = NO;
      absentCount[i] = 0;
    }
  }

  bool HasForbidden(Configuration &c, int curIter) {
    LifeState workspace;
    PutStartState(&workspace, c);

    for (int i = 0; i <= curIter + 1; i++) {
      for (int j = 0; j < numIters; j++) {
        for (int k = 0; k < forbiddenTargets[c.curs[j]].size(); k++) {
          if (Contains(&workspace, forbiddenTargets[c.curs[j]][k],
                       c.curx[j], c.cury[j]) == YES)
            return true;
        }
      }
      Run(&workspace, 1);
    }

    return false;
  }

  bool UpdateActivationCountersFail(LifeState *workspace, Configuration &conf) {
    for (int j = 0; j < numIters; j++) {
      if (Contains(workspace, conf.shiftedTargets[j]) == NO) {
        activated[j] = YES;
        absentCount[j] += MAIN_STEP;

        if (absentCount[j] > maxSurvive[conf.curs[j]]) {
          return true;
        }
      } else {
        absentCount[j] = 0;
      }
    }

    return false;
  }

  bool FilterForCurrentGenFail(LifeState *workspace) {
    for (int j = 0; j < targetFilter.size(); j++) {
      if (workspace->gen == params.filterGen[j] &&
          Contains(workspace, targetFilter[j]) == NO) {
        return true;
      }
    }

    return false;
  }

  bool IsAllActivated() {
    for (int j = 0; j < numIters; j++) {
      if (activated[j] == NO || absentCount[j] != 0) {
        return NO;
      }
    }

    return YES;
  }

  void PutStartState(LifeState *workspace, Configuration &conf) {
    ClearData(workspace);
    PutStateWSym(workspace, &conf.state, params.symmetricSearch);
    PutStateWSym(workspace, pat, params.symmetricSearch);
  }

  bool ValidateFilters(Configuration &conf) {
    LifeState workspace;
    PutStartState(&workspace, conf);

    std::vector<bool> rangeValid(params.filterGen.size(), false);

    for (int k = 0; k < params.filterGen.size(); k++)
      if (params.filterGen[k] >= 0)
        rangeValid[k] = true;

    for (int j = 0; j <= filterMaxGen; j++) {
      for (int k = 0; k < params.filterGen.size(); k++) {
        if (workspace.gen == params.filterGen[k] &&
            Contains(&workspace, targetFilter[k]) == NO)
          return false;

        if (params.filterGen[k] == -1 &&
            params.filterGenRange[k].first <= workspace.gen &&
            params.filterGenRange[k].second >= workspace.gen &&
            Contains(&workspace, targetFilter[k]) == YES)
          rangeValid[k] = true;
      }

      Run(&workspace, 1);
    }

    for (int k = 0; k < params.filterGen.size(); k++)
      if (!rangeValid[k])
        return false;

    return true;
  }

  int TestConfiguration(Configuration &conf) {
    LifeState workspace;
    ClearData(&workspace);
    PutStateWSym(&workspace, &conf.state, params.symmetricSearch);
    Join(&workspace, preIterated[conf.minIter]);
    workspace.gen = conf.minIter;

    // Initial searcher countters for absense and activation
    InitActivationCounters();

    int surviveCount = 0;

    for (int i = conf.minIter; i < iterationMaxGen; i += MAIN_STEP) {
      MainRun(workspace);

      // Fail if some catalyst is idle for too long - updates the counters for
      // them otherwise.
      if (UpdateActivationCountersFail(&workspace, conf))
        return -1;

      if (hasFilterDontReportAll) {
        // Validate filters if any of them exist. Will validate on current gen
        // of GlobalState
        if (FilterForCurrentGenFail(&workspace))
          return -1;
      }

      if (IsAllActivated()) {
        surviveCount += MAIN_STEP;
      }
      else
        surviveCount = 0;

      // If everything was actuvated and stable for stableInterval then report.
      if (surviveCount >= surviveCountForUpdate) {
        return i;
      }
    }
    return -1;
  }

  void ReportSolution(Configuration &conf, int successtime){
    LifeState workspace;
    // if reportAll - ignore filters and update fullReport
    if (reportAll) {
      ClearData(&workspace);

      PutStateWSym(&workspace, &conf.state, params.symmetricSearch);
      Copy(catalysts, &workspace);

      PutStartState(&workspace, conf);
      Copy(init, &workspace);

      Run(&workspace, successtime - surviveCountForUpdate + 2);
      Copy(afterCatalyst, &workspace);

      fullfound++;

      fullCategoryContainer->Add(init, afterCatalyst, catalysts, conf,
                                 successtime - surviveCountForUpdate + 2, 0);
    }

    // If has fitlter validate them;
    if (hasFilter) {
      if (!ValidateFilters(conf))
        return;
    }

    if (HasForbidden(conf, successtime + 3))
      return;

    // If all filters validated update results
    ClearData(&workspace);
    PutStateWSym(&workspace, &conf.state, params.symmetricSearch);
    Copy(catalysts, &workspace);

    PutStartState(&workspace, conf);
    Copy(init, &workspace);

    Run(&workspace, successtime - surviveCountForUpdate + 2);
    Copy(afterCatalyst, &workspace);

    categoryContainer->Add(init, afterCatalyst, catalysts, conf,
                           successtime - surviveCountForUpdate + 2, 0);
    found++;
}

  void SetParamsForCombine(int combineIter) {
    if (params.combineSurvive.size() - 1 < combineIter)
      surviveCountForUpdate =
          params.combineSurvive[params.combineSurvive.size() - 1];
    else
      surviveCountForUpdate = params.combineSurvive[combineIter];
  }
};

// class CategoryMultiplicator {
// public:
//   std::vector<SearchResult *> base;
//   std::vector<SearchResult *> cur;
//   CatalystSearcher *searcher;
//   int iter;

//   explicit CategoryMultiplicator(CatalystSearcher *bruteSearch) {
//     searcher = bruteSearch;
//     // searcher->categoryContainer->Sort();
//     searcher->categoryContainer->RemoveTail();

//     for (auto & category: searcher->categoryContainer->categories) {
//       base.push_back(category->results[0]);
//       cur.push_back(category->results[0]);
//     }

//     searcher->AddIterators(searcher->numIters);
//     iter = 1;
//     searcher->SetParamsForCombine(iter);
//   }

//   void CartesianMultiplication() {
//     searcher->total = base.size() * cur.size();
//     searcher->total /= 1000000;
//     searcher->idx = 0;
//     searcher->counter = 0;

//     for (auto & i : base) {
//       int idx = i->SetIters(searcher->enu, 0);
//       int baselast = i->maxGenSurvive;
//       int basefirst = i->firstGenSurvive;

//       for (auto & j : cur) {
//         searcher->IncreaseIndexAndReport(false);

//         int curlast = j->maxGenSurvive;
//         int curfirst = j->firstGenSurvive;

//         if (curlast < basefirst || baselast < curfirst)
//           continue;

//         j->SetIters(searcher->enu, idx);
//         searcher->UpdateResults();
//       }
//     }

//     // searcher->categoryContainer->Sort();
//     searcher->categoryContainer->RemoveTail();
//   }

//   void ReinitializeCurrent(int size, int iters) {
//     cur.clear();

//     for (auto & category: searcher->categoryContainer->categories) {
//       if (category
//               ->results[0]
//               ->params.size() == 3 * size)
//         cur.push_back(category->results[0]);
//     }

//     searcher->AddIterators(iters);
//     iter++;
//     searcher->SetParamsForCombine(iter);
//   }

//   void RunWithInputParams() const {
//     searcher->surviveCountForUpdate = searcher->params.stableInterval;
//     searcher->hasFilter = !searcher->params.targetFilter.empty();
//     searcher->reportAll = searcher->params.fullReportFile.length() != 0;
//     searcher->hasFilterDontReportAll =
//         searcher->hasFilter && !(searcher->reportAll);

//     CategoryContainer *found = searcher->categoryContainer;
//     searcher->categoryContainer =
//         new CategoryContainer(searcher->params.maxGen);

//     for (auto & category: found->categories) {
//       searcher->numIters =
//           category->results[0]->SetIters(searcher->enu, 0);
//       searcher->UpdateResults();
//     }
//   }
// };

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cout << "Usage CatForce.exe <in file> <nthreads>" << std::endl;
    exit(0);
  }
  int nthreads = 1;
  if (argc > 2) {
    nthreads = atoi(argv[2]);
  }
  omp_set_num_threads(nthreads);
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

  Configuration c;
  #pragma omp parallel
  {
  #pragma omp single nowait
  {
  while (!searcher.enu.done) {
    c = GetConfiguration(searcher.enu);
    // std::cout << c.curx[0] << ", " << c.cury[0] << std::endl;
    #pragma omp task default(none) shared(searcher) firstprivate(c)
    {
      int result = searcher.TestConfiguration(c);
      if(result != -1) {
        #pragma omp critical
        {
        searcher.ReportSolution(c, result);
        }
      }
      for(int i = 0; i < searcher.enu.count; i++) {
        FreeTarget(c.shiftedTargets[i]);
      }
    }
    searcher.IncreaseIndexAndReport();
    Next(searcher.enu);
  }
  }
  }
  // Print report one final time (update files with the final results).
  searcher.Report();

  // if (searcher.params.combineResults) {
  //   int startCatalysts = searcher.numIters;
  //   CategoryMultiplicator combined(&searcher);
  //   int i = 1;
  //   while (!combined.cur.empty()) {
  //     std::cout << "Combining iteration = " << i << std::endl;
  //     combined.CartesianMultiplication();

  //     std::stringstream ss;
  //     ss << "_Combined_" << ++i;

  //     searcher.Report(ss.str());
  //     combined.ReinitializeCurrent(i, startCatalysts);
  //   }

  //   std::cout << "Final result" << std::endl;

  //   combined.RunWithInputParams();

  //   searcher.Report(std::string("_Final"));
  // }

  printf("\n\nFINISH\n");
  clock_t end = clock();
  printf("Total elapsed CPU time (not wallclock if nthreads>1): %f seconds\n",
         (double)(end - searcher.begin) / CLOCKS_PER_SEC);
}
