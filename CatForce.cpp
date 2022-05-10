// CatForce - Catalyst search utility based on LifeAPI using brute force.
// Written by Michael Simkin 2015
#include "LifeAPI.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <time.h>
#include <vector>

const int none[] = {1, 0, 0, 1};
const int flipX[] = {-1, 0, 0, 1};
const int flipY[] = {1, 0, 0, -1};
const int flipXY[] = {-1, 0, 0, -1};

const int rot90clock[] = {0, 1, -1, 0};
const int rot90anti[] = {0, -1, 1, 0};
const int symmXY[] = {0, 1, 1, 0};
const int symmYX[] = {0, -1, -1, 0};

enum Symmetry {
  NONE,
  HORIZONTAL,
  HORIZONTALEVEN,
  DIAGONAL,
  DIAGONALEVENX,
  DIAGONALEVENBOTH,
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
  int searchArea[4];
  int maxW;
  int maxH;
  Symmetry symmetricSearch;
  std::vector<std::string> targetFilter;
  std::vector<int> filterdx;
  std::vector<int> filterdy;
  std::vector<int> filterGen;
  std::vector<std::pair<int, int>> filterGenRange;

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

  CatalystInput(std::string line) {
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
        forbiddenXY.push_back(std::pair<int, int>(
            atoi(elems[argi + 2].c_str()), atoi(elems[argi + 3].c_str())));

        argi += 4;
      }
    }
  }

  void Print() {
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

void ReadParams(std::string fname, std::vector<CatalystInput> &catalysts,
                SearchParams &params) {
  std::ifstream infile;
  infile.open(fname.c_str(), std::ifstream::in);

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
        catalysts.push_back(CatalystInput(line));

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
          params.filterGenRange.push_back(std::pair<int, int>(-1, -1));
        } else {
          int minGen = atoi(rangeElems[0].c_str());
          int maxGen = atoi(rangeElems[1].c_str());

          params.filterGen.push_back(-1);
          params.filterGenRange.push_back(std::pair<int, int>(minGen, maxGen));
        }

        params.targetFilter.push_back(elems[2]);
        params.filterdx.push_back(atoi(elems[3].c_str()));
        params.filterdy.push_back(atoi(elems[4].c_str()));
      }

      if (elems[0] == maxWH) {
        params.maxW = atoi(elems[1].c_str());
        params.maxH = atoi(elems[2].c_str());
      }

      if (elems[0] == combine && elems[1] == "yes") {
        params.combineResults = true;

        for (int i = 2; i < elems.size(); i++)
          params.combineSurvive.push_back(atoi(elems[i].c_str()));
      }

      if (elems[0] == symmetry) {
        if(elems[1] == "horizontal") {
          params.symmetricSearch = HORIZONTAL;
        } else if (elems[1] == "horizontaleven") {
          params.symmetricSearch = HORIZONTALEVEN;
        } else if (elems[1] == "diagonal") {
          params.symmetricSearch = DIAGONAL;
        } else if (elems[1] == "diagonalevenx") {
          params.symmetricSearch = DIAGONALEVENX;
        } else if (elems[1] == "diagonalevenboth") {
          params.symmetricSearch = DIAGONALEVENBOTH;
        }
      }


    } catch (const std::exception &ex) {
    }
  }
}

void ApplySym(LifeState *state, Symmetry sym) {
  if (sym == HORIZONTAL) {
    FlipX(state);
  } else if (sym == HORIZONTALEVEN) {
    Reverse(state->state, 0, N - 1);
  } else if (sym == DIAGONAL) {
    FlipX(state);
    FlipY(state);
  } else if (sym == DIAGONALEVENX) {
    Reverse(state->state, 0, N - 1);
    FlipY(state);
  } else if (sym == DIAGONALEVENBOTH) {
    Reverse(state->state, 0, N - 1);
    BitReverse(state);
  }
}

void PutStateWSym(LifeState *state, int x, int y, Symmetry sym) {
  PutState(state, x, y);
  if (sym == NONE)
    return;
  LifeState transformed;
  for (int i = 0; i < N; i++)
    transformed.state[i] = 0;

  transformed.min = 0;
  transformed.max = N - 1;
  transformed.gen = 0;

  Join(&transformed, state, x, y);
  ApplySym(&transformed, sym);
  PutState(&transformed);
}

void PutStateWSym(LifeState *state, Symmetry sym) {
  PutState(state);
  if (sym == NONE)
    return;

  LifeState transformed;
  Copy(&transformed, state);
  ApplySym(&transformed, sym);
  PutState(&transformed);
}

void PutStateWSym(LifeIterator *iter, Symmetry sym) {
  PutState(iter);
  if (sym == NONE)
    return;

  LifeState transformed;
  Copy(&transformed, iter->States[iter->curs], iter->curx, iter->cury);
  ApplySym(&transformed, sym);
  PutState(&transformed);
}

void GenerateStates(const std::vector<CatalystInput> &catalysts,
                    std::vector<LifeState *> &states,
                    std::vector<std::vector<LifeTarget *>> &forbidden,
                    std::vector<int> &maxSurvive) {
  for (int i = 0; i < catalysts.size(); i++) {
    std::vector<const int *> trans;
    CharToTransVec(catalysts[i].symmType, trans);

    const char *rle = catalysts[i].rle.c_str();
    int dx = catalysts[i].centerX;
    int dy = catalysts[i].centerY;
    int maxDesapear = catalysts[i].maxDesapear;

    for (int j = 0; j < trans.size(); j++) {
      int dxx = trans[j][0];
      int dxy = trans[j][1];
      int dyx = trans[j][2];
      int dyy = trans[j][3];

      states.push_back(NewState(rle, dx, dy, dxx, dxy, dyx, dyy));
      maxSurvive.push_back(maxDesapear);

      std::vector<LifeTarget *> forbidTarg;

      for (int k = 0; k < catalysts[i].forbiddenRLE.size(); k++) {
        forbidTarg.push_back(NewTarget(
            NewState(catalysts[i].forbiddenRLE[k].c_str(),
                     catalysts[i].forbiddenXY[k].first,
                     catalysts[i].forbiddenXY[k].second, dxx, dxy, dyx, dyy)));
      }

      forbidden.push_back(forbidTarg);
    }
  }
}

void InitCatalysts(std::string fname, std::vector<LifeState *> &states,
                   std::vector<std::vector<LifeTarget *>> &forbidden,
                   std::vector<int> &maxSurvive, SearchParams &params) {
  std::vector<CatalystInput> catalysts;
  ReadParams(fname, catalysts, params);
  GenerateStates(catalysts, states, forbidden, maxSurvive);
}

void XYStartGenPerState(
    const std::vector<LifeTarget *> &targets, LifeState *pat,
    const SearchParams &params, const std::vector<LifeState *> &states,
    std::vector<std::vector<std::vector<int>>> &statexyGen) {

  for (int i = 0; i < states.size(); i++) {
    std::vector<std::vector<int>> xyVec;

    for (int x = 0; x < 64; x++) {
      std::vector<int> xVec;

      for (int y = 0; y < 64; y++) {
        New();
        PutStateWSym(states[i], x, y, params.symmetricSearch);
        PutStateWSym(pat, params.symmetricSearch);
        int j;

        for (j = 0; j < params.maxGen + 5; j++) {
          if (Contains(GlobalState, targets[i], x, y) == NO) {
            break;
          }
          Run(1);
        }

        if (j == params.maxGen + 4)
          j = -1;

        xVec.push_back(j - 1);
      }

      xyVec.push_back(xVec);
    }

    statexyGen.push_back(xyVec);
  }
}

void PreIteratePat(LifeState *pat, std::vector<LifeState *> &preIterated,
                   const SearchParams &params) {
  New();
  PutStateWSym(pat, params.symmetricSearch);

  for (int i = 0; i < params.maxGen + 5; i++) {
    LifeState *t = NewState();
    Copy(t, GlobalState);
    preIterated.push_back(t);
    Run(1);
  }
}

std::string GetRLE(const std::vector<std::vector<int>> &life2d) {
  if (life2d.size() == 0)
    return "";

  if (life2d[0].size() == 0)
    return "";

  int h = life2d[0].size();

  std::stringstream result;

  int eol_count = 0;

  for (int j = 0; j < h; j++) {
    int last_val = -1;
    int run_count = 0;

    for (int i = 0; i < life2d.size(); i++) {
      int val = life2d[i][j];

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

  SearchResult(LifeState *initState, const std::vector<LifeIterator *> &iters,
               int firstGenSurviveIn, int genSurvive) {
    init = NewState();
    Copy(init, initState, COPY);

    for (int i = 0; i < iters.size(); i++) {
      params.push_back(iters[i]->curs);
      params.push_back(iters[i]->curx);
      params.push_back(iters[i]->cury);
    }

    maxGenSurvive = genSurvive;
    firstGenSurvive = firstGenSurviveIn;
  }

  int SetIters(std::vector<LifeIterator *> &iters, const int &startIdx) {
    int idx = startIdx;

    for (int i = 0; i < params.size(); i += 3) {
      iters[idx]->curs = params[i];
      iters[idx]->curx = params[i + 1];
      iters[idx]->cury = params[i + 2];
      idx++;
    }

    return idx;
  }

  void Print() {
    std::cout << "start:" << firstGenSurvive;
    std::cout << ", finish:" << maxGenSurvive << ", params: ";

    for (int i = 0; i < params.size(); i++)
      std::cout << params[i] << ",";

    std::cout << std::endl;
  }

  ~SearchResult(void) { FreeState(init); }
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
    Copy(categoryKey, catalystRemoved, COPY);
    results.push_back(firstResult);
    catDelta = catDeltaIn;
    tempCat = NewState();
    tempTest = NewState();
    maxgen = maxGen;

    ClearData(tempCat);
    Copy(tempCat, categoryKey, COPY);
    Evolve(tempCat, maxgen - tempCat->gen);
    hash = GetHash(tempCat);
  }

  void Add(SearchResult *result) { results.push_back(result); }

  bool BelongsTo(LifeState *test, const uint64_t &testHash) {
    if (testHash != hash)
      return false;

    ClearData(tempCat);
    Copy(tempCat, categoryKey, COPY);

    ClearData(tempTest);
    Copy(tempTest, test, COPY);

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

    for (int i = 0; i < results.size(); i++)
      delete results[i];
  }

  void Print() {
    for (int i = 0; i < results.size(); i++)
      results[i]->Print();
  }

  std::string RLE() {
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

    for (int l = 0; l < results.size(); l++)
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

  CategoryContainer(int maxGen) {
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
           const std::vector<LifeIterator *> &iters, int firstGenSurvive,
           int genSurvive) {
    LifeState *result = NewState();

    ClearData(result);
    Copy(result, afterCatalyst, COPY);
    Copy(result, catalysts, XOR);

    Evolve(result, maxgen - result->gen);
    uint64_t hash = GetHash(result);

    ClearData(result);
    Copy(result, afterCatalyst, COPY);
    Copy(result, catalysts, XOR);

    for (int i = 0; i < categories.size(); i++) {
      if (categories[i]->BelongsTo(result, hash)) {
        if (categories[i]->results[0]->params.size() == 3 * iters.size())
          categories[i]->Add(
              new SearchResult(init, iters, firstGenSurvive, genSurvive));
        return;
      }
    }

    LifeState *categoryKey = NewState();
    Copy(categoryKey, afterCatalyst, COPY);
    Copy(categoryKey, catalysts, XOR);

    categories.push_back(new Category(
        categoryKey, new SearchResult(init, iters, firstGenSurvive, genSurvive),
        catDelta, maxgen));
  }

  void Sort() {
    for (int i = 0; i < categories.size(); i++)
      categories[i]->Sort();
  }

  void Print() {
    for (int i = 0; i < categories.size(); i++)
      categories[i]->Print();
  }

  void RemoveTail() {
    for (int i = 0; i < categories.size(); i++)
      categories[i]->RemoveTail();
  }

  std::string CategoriesRLE() {
    std::stringstream ss;
    for (int i = 0; i < categories.size(); i++) {
      ss << categories[i]->RLE();
    }

    return ss.str();
  }
};

class CatalystSearcher {
public:
  std::string result;
  clock_t begin;
  std::vector<LifeState *> states;
  std::vector<int> maxSurvive;
  SearchParams params;
  LifeState *pat;
  int numIters;
  std::vector<LifeIterator *> iters;
  std::vector<LifeTarget *> targetFilter;
  std::vector<LifeTarget *> targets;
  std::vector<LifeTarget *> preShifted;
  std::vector<std::vector<LifeTarget *>> forbiddenTargets;
  std::vector<std::vector<std::vector<int>>> statexyGen;
  std::vector<LifeState *> preIterated;
  std::vector<int> activated;
  std::vector<int> absentCount;
  std::vector<LifeState *> statesArr;
  clock_t current;
  long long idx;
  int found;
  long long total;
  std::string fullReport;
  unsigned short int counter;
  CategoryContainer *categoryContainer;

  // flags and memeber for the search

  bool hasFilter;
  bool reportAll;
  bool hasFilterDontReportAll;

  int filterMaxGen;
  int iterationMaxGen;

  int surviveCountForUpdate;

  LifeState *init;
  LifeState *afterCatalyst;
  LifeState *catalysts;

  void Init(char *inputFile) {
    result = "x = 0, y = 0, rule = B3/S23\n";
    begin = clock();
    InitCatalysts(inputFile, states, forbiddenTargets, maxSurvive, params);
    pat = NewState(params.pat.c_str(), params.xPat, params.yPat);
    numIters = params.numCatalysts;
    categoryContainer = new CategoryContainer(params.maxGen);

    for (int i = 0; i < params.targetFilter.size(); i++)
      targetFilter.push_back(NewTarget(params.targetFilter[i].c_str(),
                                       params.filterdx[i], params.filterdy[i]));

    for (int i = 0; i < states.size(); i++)
      targets.push_back(NewTarget(states[i]));

    XYStartGenPerState(targets, pat, params, states, statexyGen);

    PreIteratePat(pat, preIterated, params);
    AddIterators(numIters);

    current = clock();
    idx = 0;
    found = 0;
    total = 1;
    counter = 0;

    int fact = 1;

    for (int i = 0; i < numIters; i++) {
      total *= (iters[i]->w);
      total *= (iters[i]->h);
      total *= (iters[i]->s);
      fact *= (i + 1);
    }

    total /= fact;

    std::cout << "Approximated Total: " << total << std::endl;
    total = total / 1000000;

    if (total == 0)
      total++;

    fullReport = "x = 0, y = 0, rule = B3/S23\n";

    hasFilter = params.targetFilter.size() > 0;
    reportAll = params.fullReportFile.length() != 0;
    hasFilterDontReportAll = hasFilter && !reportAll;

    filterMaxGen = FilterMaxGen();
    iterationMaxGen = params.maxGen;
    numIters = numIters;

    init = NewState();
    afterCatalyst = NewState();
    catalysts = NewState();

    surviveCountForUpdate = params.stableInterval;

    if (params.combineResults) {
      if (params.combineSurvive.size() == 0)
        params.combineSurvive.push_back(1);

      surviveCountForUpdate = params.combineSurvive[0];
      hasFilter = false;
      hasFilterDontReportAll = false;
    }
  }

  void AddIterators(int num) {
    for (int i = 0; i < num; i++) {
      iters.push_back(NewIterator(&states[0], params.searchArea[0],
                                  params.searchArea[1], params.searchArea[2],
                                  params.searchArea[3], states.size()));
      activated.push_back(0);
      absentCount.push_back(0);
      preShifted.push_back(NewTarget(NewState(), NewState()));
    }

    numIters = iters.size();
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

  void Report(std::string suffix) {
    std::string temp = params.outputFile;
    params.outputFile = params.outputFile + suffix + std::string(".rle");
    Report();
    params.outputFile = temp;
  }

  void Report(bool saveFile = true) {
    float percent = (idx / 10000) / (total * 1.0);
    int sec = (clock() - begin) / CLOCKS_PER_SEC + 1;
    int estimation = 0;
    int checkPerSecond = idx / (sec * 1000);

    if (percent > 0)
      estimation = (sec * 100) / percent;

    std::cout << std::setprecision(1) << std::fixed << percent << "%,"
              << idx / 1000000 << "M/" << total
              << "M, cats/find: " << categoryContainer->categories.size() << "/"
              << found << ", now: ";
    PrintTime(sec);
    std::cout << ", est: ";
    PrintTime(estimation);
    std::cout << ", " << std::setprecision(1) << std::fixed << checkPerSecond
              << "K/sec" << std::endl;

    // Linear report seems to be non needed for now. Don't delete might have
    // some use in the future. std::ofstream
    // resultsFile(params.outputFile.c_str()); resultsFile << result;
    // resultsFile.close();

    categoryContainer->Sort();

    if (saveFile) {
      std::cout << "Saving " << params.outputFile << std::endl;

      std::ofstream catResultsFile(params.outputFile.c_str());
      catResultsFile << "x = 0, y = 0, rule = B3/S23\n";
      catResultsFile << categoryContainer->CategoriesRLE();
      catResultsFile.close();

      if (params.fullReportFile.length() != 0) {
        std::ofstream allfile(params.fullReportFile.c_str());
        allfile << fullReport;
        allfile.close();
      }
    }
  }

  void PrintTime(int sec) {
    int hr = sec / 3600;
    int min = (sec / 60) - hr * 60;
    int secs = sec - 3600 * hr - 60 * min;
    std::cout << hr << ":" << min << ":" << secs;
    return;
  }
  void IncreaseIndexAndReport(bool saveFile = true) {
    counter++;

    if (counter == 0) {
      idx += 65536;

      if (idx % (1048576) == 0) {
        if ((double)(clock() - current) / CLOCKS_PER_SEC > 10) {
          current = clock();
          Report(saveFile);
        }
      }
    }
  }

  int ValidateMinWidthHeight() {
    // Fast path: we know that the iterators are ordered by x.
    if(iters[0]->curx - iters[numIters-1]->curx >= params.maxW)
      return NO;

    int minY = iters[0]->cury;
    int maxY = iters[0]->cury;

    for (int i = 1; i < numIters; i++) {
      if (iters[i]->cury > maxY)
        maxY = iters[i]->cury;

      if (iters[i]->cury < minY)
        minY = iters[i]->cury;
    }

    if (maxY - minY >= params.maxH)
      return NO;
    else
      return YES;
  }

  int LastNonActiveGeneration() {
    int minIter = statexyGen[iters[0]->curs][(iters[0]->curx + 64) % 64]
                            [(iters[0]->cury + 64) % 64];

    for (int i = 1; i < numIters; i++) {
      int startGen = statexyGen[iters[i]->curs][(iters[i]->curx + 64) % 64]
                               [(iters[i]->cury + 64) % 64];

      if (startGen < minIter)
        minIter = startGen;

      if (minIter < params.startGen)
        break;
    }

    return minIter;
  }

  void PreShiftIters() {
    for (int i = 0; i < numIters; i++) {
      Copy(preShifted[i]->wanted,   targets[iters[i]->curs]->wanted, iters[i]->curx, iters[i]->cury);
      Copy(preShifted[i]->unwanted, targets[iters[i]->curs]->unwanted, iters[i]->curx, iters[i]->cury);
    }
  }

  void PutItersState() {
    for (int i = 0; i < numIters; i++) {
      PutStateWSym(iters[i], params.symmetricSearch);
    }
  }

  int CatalystCollide() {
    Run(1);

    for (int i = 0; i < numIters; i++) {
      if (Contains(GlobalState, preShifted[i]) == NO) {
        return YES;
      }
    }

    return NO;
  }

  void InitActivationCounters() {
    for (int i = 0; i < numIters; i++) {
      activated[i] = NO;
      absentCount[i] = 0;
    }
  }

  bool HasForbidden(int curIter) {
    PutCurrentState();

    for (int i = 0; i <= curIter + 1; i++) {
      for (int j = 0; j < numIters; j++) {
        for (int k = 0; k < forbiddenTargets[iters[j]->curs].size(); k++) {
          if (Contains(GlobalState, forbiddenTargets[iters[j]->curs][k],
                       iters[j]->curx, iters[j]->cury) == YES)
            return true;
        }
      }
      Run(1);
    }

    return false;
  }

  bool UpdateActivationCountersFail() {
    for (int j = 0; j < numIters; j++) {
      if (Contains(GlobalState, preShifted[j]) == NO) {
        activated[j] = YES;
        absentCount[j]++;

        if (absentCount[j] > maxSurvive[iters[j]->curs]) {
          return true;
        }
      } else {
        absentCount[j] = 0;
      }
    }

    return false;
  }

  bool FilterForCurrentGenFail() {
    for (int j = 0; j < targetFilter.size(); j++) {
      if (GlobalState->gen == params.filterGen[j] &&
          Contains(GlobalState, targetFilter[j]) == NO) {
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

  void PutCurrentState() {
    New();
    PutItersState();
    PutStateWSym(pat, params.symmetricSearch);
  }

  bool ValidateFilters(const int &maxFilterGen) {
    std::vector<bool> rangeValid(params.filterGen.size(), false);

    for (int k = 0; k < params.filterGen.size(); k++)
      if (params.filterGen[k] >= 0)
        rangeValid[k] = true;

    for (int j = 0; j <= maxFilterGen; j++) {
      for (int k = 0; k < params.filterGen.size(); k++) {
        if (GlobalState->gen == params.filterGen[k] &&
            Contains(GlobalState, targetFilter[k]) == NO)
          return false;

        if (params.filterGen[k] == -1 &&
            params.filterGenRange[k].first <= GlobalState->gen &&
            params.filterGenRange[k].second >= GlobalState->gen &&
            Contains(GlobalState, targetFilter[k]) == YES)
          rangeValid[k] = true;
      }

      Run(1);
    }

    for (int k = 0; k < params.filterGen.size(); k++)
      if (!rangeValid[k])
        return false;

    return true;
  }

  void UpdateResults() {

    int minIter = LastNonActiveGeneration();

    // Activation before first generation allowed to be activated
    if (minIter < params.startGen || minIter >= params.lastGen)
      return;

    // Place catalysts first and check if they collide.
    New();

    PreShiftIters();
    PutItersState();
    if (CatalystCollide() == YES)
      return;

    PutState(preIterated[minIter]);

    // Initial searcher countters for absense and activation
    InitActivationCounters();

    int surviveCount = 0;

    for (int i = minIter; i < iterationMaxGen; i++) {
      Run(1);

      // Fail if some catalyst is idle for too long - updates the counters for
      // them otherwise.
      if (UpdateActivationCountersFail())
        break;

      // const bool optimization - will skip this always.
      if (hasFilterDontReportAll) {
        // Validate filters if any of them exist. Will validate on current gen
        // of GlobalState
        if (FilterForCurrentGenFail())
          break;
      }

      if (IsAllActivated())
        surviveCount++;
      else
        surviveCount = 0;

      // If everything was actuvated and stable for stableInterval then report.
      if (surviveCount >= surviveCountForUpdate) {
        bool valid = true;

        int genSurvive;

        for (genSurvive = i; genSurvive < iterationMaxGen; genSurvive++) {
          Run(1);

          if (UpdateActivationCountersFail())
            break;
        }
        // If has fitlter validate them;
        if (hasFilter) {
          PutCurrentState();
          valid = ValidateFilters(filterMaxGen);
        }

        if (!valid)
          break;

        if (HasForbidden(i + 3))
          break;

        // If all filters validated update results
        if (valid) {
          ClearData(catalysts);
          ClearData(init);
          ClearData(afterCatalyst);

          New();
          PutItersState();

          Copy(catalysts, GlobalState, COPY);
          PutCurrentState();
          Copy(init, GlobalState, COPY);
          Run(i - surviveCountForUpdate + 2);
          Copy(afterCatalyst, GlobalState, COPY);

          categoryContainer->Add(init, afterCatalyst, catalysts, iters,
                                 i - surviveCountForUpdate + 2, genSurvive);
          PutCurrentState();
          result.append(GetRLE(GlobalState));
          result.append("100$");
          found++;
        }

        // if reportAll - ignore filters and update fullReport
        if (reportAll) {
          PutCurrentState();

          fullReport.append(GetRLE(GlobalState));
          fullReport.append("100$");
        }
        break;
      }
    }
  }

  void SetParamsForCombine(int combineIter) {
    if (params.combineSurvive.size() - 1 < combineIter)
      surviveCountForUpdate =
          params.combineSurvive[params.combineSurvive.size() - 1];
    else
      surviveCountForUpdate = params.combineSurvive[combineIter];
  }
};

class CategoryMultiplicator {
public:
  std::vector<SearchResult *> base;
  std::vector<SearchResult *> cur;
  CatalystSearcher *searcher;
  int iter;

  CategoryMultiplicator(CatalystSearcher *bruteSearch) {
    searcher = bruteSearch;
    searcher->categoryContainer->Sort();
    searcher->categoryContainer->RemoveTail();

    for (int i = 0; i < searcher->categoryContainer->categories.size(); i++) {
      base.push_back(searcher->categoryContainer->categories[i]->results[0]);
      cur.push_back(searcher->categoryContainer->categories[i]->results[0]);
    }

    searcher->AddIterators(searcher->numIters);
    iter = 1;
    searcher->SetParamsForCombine(iter);
  }

  void CartesianMultiplication() {
    searcher->total = base.size() * cur.size();
    searcher->total /= 1000000;
    searcher->idx = 0;
    searcher->counter = 0;

    for (int i = 0; i < base.size(); i++) {
      int idx = base[i]->SetIters(searcher->iters, 0);
      int baselast = base[i]->maxGenSurvive;
      int basefirst = base[i]->firstGenSurvive;

      for (int j = 0; j < cur.size(); j++) {
        searcher->IncreaseIndexAndReport(false);

        int curlast = cur[j]->maxGenSurvive;
        int curfirst = cur[j]->firstGenSurvive;

        if (curlast < basefirst || baselast < curfirst)
          continue;

        cur[j]->SetIters(searcher->iters, idx);
        searcher->UpdateResults();
      }
    }

    searcher->categoryContainer->Sort();
    searcher->categoryContainer->RemoveTail();
  }

  void ReinitializeCurrent(int size, int iters) {
    cur.clear();

    for (int i = 0; i < searcher->categoryContainer->categories.size(); i++) {
      if (searcher->categoryContainer->categories[i]
              ->results[0]
              ->params.size() == 3 * size)
        cur.push_back(searcher->categoryContainer->categories[i]->results[0]);
    }

    searcher->AddIterators(iters);
    iter++;
    searcher->SetParamsForCombine(iter);
  }

  void RunWithInputParams() {
    searcher->surviveCountForUpdate = searcher->params.stableInterval;
    searcher->hasFilter = searcher->params.targetFilter.size() > 0;
    searcher->reportAll = searcher->params.fullReportFile.length() != 0;
    searcher->hasFilterDontReportAll =
        searcher->hasFilter && !(searcher->reportAll);

    CategoryContainer *found = searcher->categoryContainer;
    searcher->categoryContainer =
        new CategoryContainer(searcher->params.maxGen);

    for (int i = 0; i < found->categories.size(); i++) {
      searcher->numIters =
          found->categories[i]->results[0]->SetIters(searcher->iters, 0);
      searcher->UpdateResults();
    }
  }
};

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cout << "Usage CatForce.exe <in file>" << std::endl;
    exit(0);
  }

  // LifeAPI initialization.
  New();

  std::cout << "Input: " << argv[1] << std::endl
            << "Initializing please wait..." << std::endl;

  CatalystSearcher searcher;
  searcher.Init(argv[1]);

  std::cout << std::endl
            << "Initialization finished, searching..." << std::endl
            << std::endl;

  const bool validateWH = searcher.params.maxW > 0 && searcher.params.maxH > 0;
  const int numIters = searcher.numIters;
  // Main loop of search on iters
  do {
    int valid = YES;

    // width-height validation enabled
    if (validateWH) {
      valid = searcher.ValidateMinWidthHeight();
    }

    searcher.IncreaseIndexAndReport();

    // Valid remains YES after width-height Validation
    if (valid == NO)
      continue;

    searcher.UpdateResults();

  } while (Next(&searcher.iters[0], numIters, NO));

  // Print report one final time (update files with the final results).
  searcher.Report();

  if (searcher.params.combineResults) {
    int startCatalysts = searcher.numIters;
    CategoryMultiplicator combined(&searcher);
    int i = 1;
    while (combined.cur.size() > 0) {
      std::cout << "Combining iteration = " << i << std::endl;
      combined.CartesianMultiplication();

      std::stringstream ss;
      ss << "_Combined_" << ++i;

      searcher.Report(ss.str());
      combined.ReinitializeCurrent(i, startCatalysts);
    }

    std::cout << "Final result" << std::endl;

    combined.RunWithInputParams();

    searcher.Report(std::string("_Final"));
  }

  printf("\n\nFINISH\n");
  clock_t end = clock();
  printf("Total elapsed time: %f seconds\n",
         (double)(end - searcher.begin) / CLOCKS_PER_SEC);
}
