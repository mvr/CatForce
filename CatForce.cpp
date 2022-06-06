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
  int maxDesapear;
  int centerX;
  int centerY;
  char symmType;
  std::vector<std::string> forbiddenRLE;
  std::vector<std::pair<int, int>> forbiddenXY;
  std::string requiredRLE;
  std::pair<int, int> requiredXY;

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
    maxDesapear = atoi(elems[2].c_str());
    centerX = atoi(elems[3].c_str());
    centerY = atoi(elems[4].c_str());
    symmType = elems[5].at(0);

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
      } else {
        std::cout << "Unknown catalyst attribute: " << elems[argi] << std::endl;
        exit(1);
      }
    }
  }

  void Print() const {
    std::cout << rle << " " << maxDesapear << " " << centerX << " " << centerY
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
          } else if ( rest == "/" or rest == "/odd") {
            params.symmetryChain = {ReflectYeqNegX};
          } else if ( rest == "\\" or rest == "\\odd") {
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

struct Configuration {
  int count;
  std::array<int, MAX_CATALYSTS> curx;
  std::array<int, MAX_CATALYSTS> cury;
  std::array<int, MAX_CATALYSTS> curs;
  // int minIter;
  LifeState state;
  LifeState catalystsState;
  // std::vector<LifeTarget> shiftedTargets;
};

// struct Enumerator {
//   int done;
//   const int count;

//   const int x;
//   const int y;
//   const int w;
//   const int h;
//   const int s;

//   const int maxW;
//   const int maxH;

//   const int startGen;
//   const int lastGen;
//   const int maxGen;

//   const std::vector<LifeState> states;
//   const std::vector<LifeTarget> targets;
//   const std::vector<std::vector<std::vector<int> > >  activations;

//   std::vector<int> curx;
//   std::vector<int> cury;
//   std::vector<int> curs;
//   std::vector<int> cumulMinX;
//   std::vector<int> cumulMaxX;
//   std::vector<int> cumulMinY;
//   std::vector<int> cumulMaxY;
//   std::vector<LifeTarget> shiftedTargets;
//   std::vector<LifeState> cumulative; // backwards! e.g. 1|2|3, 2|3, 3
//   std::vector<int> cumulActivation;

//   Enumerator(SearchParams &params, std::vector<LifeState> &states, std::vector<LifeTarget> &targets, std::vector<std::vector<std::vector<int> > > &activations) :
//       count(params.numCatalysts),
//       x(params.searchArea[0]),
//       y(params.searchArea[1]),
//       w(params.searchArea[2]),
//       h(params.searchArea[3]),
//       s(states.size()),
//       maxW(params.maxW),
//       maxH(params.maxH),
//       startGen(params.startGen),
//       lastGen(params.lastGen),
//       maxGen(params.maxGen),
//       states(states),
//       targets(targets),
//       activations(activations)
//   {
//     done = false;
//     curx = std::vector<int>(count, x);
//     cury = std::vector<int>(count, y);
//     curs = std::vector<int>(count, 0);
//     cumulMinX = std::vector<int>(count, x);
//     cumulMaxX = std::vector<int>(count, x);
//     cumulMinY = std::vector<int>(count, y);
//     cumulMaxY = std::vector<int>(count, y);
//     cumulActivation = std::vector<int>(count, activations[0][(x + 64) % 64][(y + 64) % 64]);

//     for (int i = 0; i < count; i++) {
//       shiftedTargets.push_back(LifeTarget());
//       cumulative.push_back(LifeState());
//     }
//     for (int i = count - 1; i >= 0; i--) {
//       shiftedTargets[i].wanted.Copy(targets[curs[i]].wanted, curx[i], cury[i]);
//       shiftedTargets[i].unwanted.Copy(targets[curs[i]].unwanted, curx[i], cury[i]);
//       cumulative[i].Copy(shiftedTargets[i].wanted);
//       if(i < count-1)
//         cumulative[i].Join(cumulative[i + 1]);
//     }
//   }

//   int Next(const int i) {
//     while (true) {
//       int n = NaiveNext(i);
//       if(n == FAIL) {
//         if(i == 0)
//           done = true;
//         return FAIL;
//       }

//       const int thisactivation = activations[curs[i]][(curx[i] + 64) % 64][(cury[i] + 64) % 64];

//       if (i == count - 1) {
//         // Check collision time (only relevant condition)
//         if (thisactivation < startGen || thisactivation > lastGen){
//           continue;
//         }

//         shiftedTargets[i].wanted.Copy(targets[curs[i]].wanted, curx[i], cury[i]);

//         break;
//       }

//       const int lastactivation = cumulActivation[i+1];

//       if(thisactivation < lastactivation)
//         continue;

//       if(thisactivation == lastactivation) { // Then break the tie
//         if (curx[i] < curx[i + 1])
//           continue;
//         if (curx[i] == curx[i + 1]) {
//           if (cury[i] < cury[i + 1])
//             continue;
//           if (cury[i] == cury[i + 1]) {
//             if (curs[i] <= curs[i + 1]) {
//               continue;
//             }
//           }
//         }
//       }

//       // Check bounds
//       if (maxW != -1)
//         if( curx[i] - cumulMinX[i+1] > maxW ||
//             cumulMaxX[i+1] - curx[i] > maxW ) {
//           continue;
//         }

//       if (maxH != -1) {
//         if( cury[i] - cumulMinY[i+1] > maxH ||
//             cumulMaxY[i+1] - cury[i] > maxH ) {
//           continue;
//         }
//       }

//       shiftedTargets[i].wanted.Copy(targets[curs[i]].wanted, curx[i], cury[i]);

//       // Check overlap
//       LifeState temp = cumulative[i + 1];
//       temp.Join(shiftedTargets[i].wanted);
//       temp.Step();
//       if (!temp.Contains(shiftedTargets[i].wanted))
//         continue;

//       break;
//     }
//     // Not needed in the loop
//     shiftedTargets[i].unwanted.Copy(targets[curs[i]].unwanted, curx[i], cury[i]);

//     cumulActivation[i] = activations[curs[i]][(curx[i] + 64) % 64][(cury[i] + 64) % 64];

//     if(i == count-1) {
//       cumulMinX[i] = curx[i];
//       cumulMaxX[i] = curx[i];
//       cumulMinY[i] = cury[i];
//       cumulMaxY[i] = cury[i];

//       cumulative[i].Copy(shiftedTargets[i].wanted);
//     }


//     if(i < count-1) {
//       // Update bounds
//       cumulMinX[i] = std::min(cumulMinX[i + 1], curx[i]);
//       cumulMaxX[i] = std::max(cumulMaxX[i + 1], curx[i]);
//       cumulMinY[i] = std::min(cumulMinY[i + 1], cury[i]);
//       cumulMaxY[i] = std::max(cumulMaxY[i + 1], cury[i]);

//       // Update cumulative
//       cumulative[i].Copy(shiftedTargets[i].wanted);
//       cumulative[i].Join(cumulative[i + 1]);
//     }

//     return SUCCESS;
//   }

//   int NaiveNext(const int i) {
//     if(i == count)
//       return FAIL;

//     curs[i]++;
//     if (curs[i] < s)
//       return SUCCESS;
//     curs[i] = 0;

//     cury[i]++;
//     if (cury[i] < y + h)
//       return SUCCESS;
//     cury[i] = y;

//     curx[i]++;
//     if (curx[i] < x + w)
//       return SUCCESS;
//     curx[i] = x;

//     return Next(i+1);
//   }

//   int Next() {
//     return Next(0);
//   }

//   Configuration GetConfiguration() {
//     Configuration c;
//     c.count = count;
//     // c.curx = curx;
//     // c.cury = cury;
//     // c.curs = curs;
//     // c.minIter = cumulActivation[count-1];
//     c.state = cumulative[0];
//     c.shiftedTargets = shiftedTargets;
//     // for(int i = 0; i < enu.count; i++) {
//     //   c.shiftedTargets.push_back(NewTarget(enu.shiftedTargets[i].wanted, enu.shiftedTargets[i].unwanted));
//     // }
//     return c;
//   }
// };

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

void GenerateStates(const std::vector<CatalystInput> &catalysts,
                    std::vector<LifeState> &states,
                    std::vector<LifeState> &required,
                    std::vector<std::vector<LifeTarget>> &forbidden,
                    std::vector<int> &maxMissing) {

  for (const auto & catalyst : catalysts) {
    std::vector<SymmetryTransform> trans;
    CharToTransVec(catalyst.symmType, trans);

    const char *rle = catalyst.rle.c_str();
    int dx = catalyst.centerX;
    int dy = catalyst.centerY;
    int maxDesapear = catalyst.maxDesapear;

    for (auto & tran : trans) {
      states.push_back(LifeState::Parse(rle, dx, dy, tran));
      maxMissing.push_back(maxDesapear);

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
    }
  }
}

void InitCatalysts(const std::string &fname, std::vector<LifeState> &states,
                   std::vector<LifeState> &required,
                   std::vector<std::vector<LifeTarget>> &forbidden,
                   std::vector<int> &maxMissing, SearchParams &params) {
  std::vector<CatalystInput> catalysts;
  ReadParams(fname, catalysts, params);
  GenerateStates(catalysts, states, required, forbidden, maxMissing);
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
  LifeState init;

  // iters state in form of integers
  // std::vector<int> params;
  int maxGenSurvive;
  int firstGenSurvive;

  SearchResult(LifeState &initState, const Configuration &conf,
               int firstGenSurviveIn, int genSurvive) {
    init.Copy(initState);

    // for (int i = 0; i < conf.curs.size(); i++) {
    //   params.push_back(conf.curs[i]);
    //   params.push_back(conf.curx[i]);
    //   params.push_back(conf.cury[i]);
    // }

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
  std::vector<LifeState> catalysts;
  std::vector<int> maxMissing;
  SearchParams params;
  LifeState pat;
  int numIters{};
  // Enumerator *enu;
  std::vector<LifeTarget> targetFilter;
  std::vector<LifeTarget> targets;
  std::vector<std::vector<LifeTarget>> forbiddenTargets;
  std::vector<LifeState> requiredParts;
  std::vector<std::vector<LifeState>> catalystCollisionMasks;
  std::vector<LifeState> catalystReactionMasks;
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

  void Init(const char *inputFile, int nthreads) {
    begin = clock();
    InitCatalysts(inputFile, catalysts, requiredParts, forbiddenTargets, maxMissing, params);
    pat = LifeState::Parse(params.pat.c_str(), params.xPat, params.yPat);
    numIters = params.numCatalysts;
    categoryContainer = new CategoryContainer(params.maxGen);
    fullCategoryContainer = new CategoryContainer(params.maxGen);

    for (auto & state : catalysts)
      targets.push_back(LifeTarget(state));

    for (int i = 0; i < params.targetFilter.size(); i++)
      targetFilter.push_back(LifeTarget::Parse(params.targetFilter[i].c_str(),
                                               params.filterdx[i], params.filterdy[i]));

    catalystCollisionMasks = std::vector<std::vector<LifeState>>(
        catalysts.size(), std::vector<LifeState>(catalysts.size()));
    catalystReactionMasks = std::vector<LifeState>(catalysts.size());

    for (int s = 0; s < catalysts.size(); s++) {
      catalystReactionMasks[s] = catalysts[s].BigZOI();
      catalystReactionMasks[s].Transform(Rotate180OddBoth);
      for (int t = 0; t < catalysts.size(); t++) {
        catalystCollisionMasks[s][t] = CollisionMask(catalysts[s], catalysts[t]);
      }
    }

    // XYStartGenPerState(targets, pat, params, catalysts, statexyGen, nthreads);

    // enu = new Enumerator(params, catalysts, targets, statexyGen);

    // PreIteratePat(pat, preIterated, params);
    // AddIterators(numIters);

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

  bool HasForbidden(Configuration &c, int curIter) {
    LifeState workspace;
    PutStartState(workspace, c);

    for (int i = 0; i <= curIter + 1; i++) {
      for (int j = 0; j < numIters; j++) {
        for (int k = 0; k < forbiddenTargets[c.curs[j]].size(); k++) {
          if (workspace.Contains(forbiddenTargets[c.curs[j]][k],
                       c.curx[j], c.cury[j]) == true)
            return true;
        }
      }
      workspace.Step();
    }

    return false;
  }

  // bool UpdateActivationCountersFail(LifeState &workspace, Configuration &conf) {
  //   for (int j = 0; j < numIters; j++) {
  //     if (workspace.Contains(conf.shiftedTargets[j]) == false) {
  //       activated[j] = true;
  //       absentCount[j] += MAIN_STEP;

  //       if (absentCount[j] > maxMissing[conf.curs[j]]) {
  //         return true;
  //       }
  //     } else {
  //       absentCount[j] = 0;
  //     }
  //   }

  //   return false;
  // }

  bool FilterForCurrentGenFail(LifeState &workspace) {
    for (int j = 0; j < targetFilter.size(); j++) {
      if (workspace.gen == params.filterGen[j] &&
          workspace.Contains(targetFilter[j]) == false) {
        return true;
      }
    }

    return false;
  }

  void PutStartState(LifeState &workspace, Configuration &conf) {
    workspace.Clear();
    workspace.JoinWSymChain(conf.state, params.symmetryChain);
    workspace.JoinWSymChain(pat, params.symmetryChain);
  }

  bool ValidateFilters(Configuration &conf) {
    LifeState workspace;
    PutStartState(workspace, conf);

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

  // int TestConfiguration(Configuration &conf) {
  //   if(conf.minIter == -1)
  //     return -1; // Temporary fix

  //   LifeState workspace;
  //   workspace.JoinWSymChain(conf.state, params.symmetryChain);
  //   workspace.Join(preIterated[conf.minIter]);
  //   workspace.gen = conf.minIter;

  //   // Initial searcher countters for absense and activation
  //   InitActivationCounters();

  //   int surviveCount = 0;

  //   for (int i = conf.minIter; i < params.maxGen; i += MAIN_STEP) {
  //     MainRun(workspace);

  //     // Fail if some catalyst is idle for too long - updates the counters for
  //     // them otherwise.
  //     if (UpdateActivationCountersFail(workspace, conf))
  //       return -1;

  //     if (hasFilter && !reportAll) {
  //       // Validate filters if any of them exist. Will validate on current gen
  //       // of GlobalState
  //       if (FilterForCurrentGenFail(workspace))
  //         return -1;
  //     }

  //     if (IsAllActivated()) {
  //       surviveCount += MAIN_STEP;
  //     }
  //     else
  //       surviveCount = 0;

  //     // If everything was actuvated and stable for stableInterval then report.
  //     if (surviveCount >= params.stableInterval) {
  //       return i;
  //     }
  //   }
  //   return -1;
  // }

  void ReportSolution(Configuration &conf, int successtime){
    LifeState init;
    LifeState afterCatalyst;
    LifeState catalysts;

    LifeState workspace;
    // if reportAll - ignore filters and update fullReport
    if (reportAll) {
      workspace.JoinWSymChain(conf.state, params.symmetryChain);
      catalysts.Copy(workspace);

      PutStartState(workspace, conf);
      init.Copy(workspace);

      workspace.Step(successtime - params.stableInterval + 2);
      afterCatalyst.Copy(workspace);

      fullfound++;

      fullCategoryContainer->Add(init, afterCatalyst, catalysts, conf,
                                 successtime - params.stableInterval + 2, 0);
    }

    // If has fitlter validate them;
    if (hasFilter) {
      if (!ValidateFilters(conf))
        return;
    }

    if (HasForbidden(conf, successtime + 3))
      return;

    // If all filters validated update results
    workspace.Clear();
    workspace.JoinWSymChain(pat, params.symmetryChain);
    workspace.Join(conf.catalystsState);
    init.Copy(workspace);

    workspace.Step(successtime - params.stableInterval + 2);
    afterCatalyst.Copy(workspace);

    categoryContainer->Add(init, afterCatalyst, conf.catalystsState, conf,
                           successtime - params.stableInterval + 2, 0);
    found++;
  }

  // void Search() {
  //   Configuration c;
  //   while (!enu->done) {
  //     c = enu->GetConfiguration();
  //     int result = TestConfiguration(c);
  //     if(result != -1) {
  //       ReportSolution(c, result);
  //     }
  //     IncreaseIndexAndReport();
  //     enu->Next();
  //   }
  // }

  void Search() {
    Configuration config;
    config.count = 0;
    config.state.JoinWSymChain(pat, params.symmetryChain);
    LifeState history = config.state;

    LifeState bounds =
        LifeState::SolidRect(params.searchArea[0], params.searchArea[1],
                             params.searchArea[2], params.searchArea[3]);

    std::vector<LifeState> masks(catalysts.size());
    for (int s = 0; s < catalysts.size(); s++) {
      masks[s] = pat.Convolve(catalystReactionMasks[s]);
      masks[s].Copy(bounds, ORNOT);
    }

    std::vector<LifeTarget> shiftedTargets(params.numCatalysts);

    RecursiveSearch(config, history, LifeState(), masks, shiftedTargets,
                    std::array<int, MAX_CATALYSTS>(), std::array<int, MAX_CATALYSTS>(),
                    std::array<bool, MAX_CATALYSTS>(), std::array<bool, MAX_CATALYSTS>());
  }

  void
  RecursiveSearch(Configuration config, LifeState history, LifeState required,
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

        if (missingTime[i] > maxMissing[i])
          return;
      }


      // Try adding a catalyst
      if (config.state.gen >= params.startGen && config.count != params.numCatalysts) {
        LifeState newcells = config.state;
        newcells.Copy(history, ANDNOT);
        history.Copy(config.state, OR);

        if (!newcells.IsEmpty()) {
          for(int s = 0; s < catalysts.size(); s++) {
            LifeState newPlacements = newcells.Convolve(catalystReactionMasks[s]);
            newPlacements.Copy(masks[s], ANDNOT);

            while (!newPlacements.IsEmpty()) {
              // Do the placement
              auto newPlacement = newPlacements.FirstOn();

              if (config.count == 0) {
                std::cout << "Placing catalyst " << s << " at "
                          << newPlacement.first << ", " << newPlacement.second
                          << std::endl;
                // newPlacements.Print();
              }

              Configuration newConfig = config;
              newConfig.count += 1;
              newConfig.curx[config.count] = newPlacement.first;
              newConfig.cury[config.count] = newPlacement.second;
              newConfig.curs[config.count] = s;

              LifeState shiftedCatalyst = catalysts[s];
              shiftedCatalyst.Move(newPlacement.first, newPlacement.second);
              newConfig.state.JoinWSymChain(shiftedCatalyst, params.symmetryChain);
              newConfig.catalystsState.JoinWSymChain(shiftedCatalyst, params.symmetryChain);

              shiftedTargets[config.count] = LifeTarget(shiftedCatalyst);

              LifeState newRequired = required;
              newRequired.Join(requiredParts[s], newPlacement.first, newPlacement.second);

              std::array<int, MAX_CATALYSTS> newMissingTime = missingTime;
              std::array<int, MAX_CATALYSTS> newRecoveredTime = recoveredTime;
              std::array<bool, MAX_CATALYSTS> newHasReacted = hasReacted;
              std::array<bool, MAX_CATALYSTS> newHasRecovered = hasRecovered;

              std::vector<LifeState> newMasks = masks;

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

              RecursiveSearch(newConfig, history, newRequired, newMasks, shiftedTargets, newMissingTime,
                              newRecoveredTime, newHasReacted, newHasRecovered);

              masks[s].Set(newPlacement.first, newPlacement.second);
              newPlacements.Erase(newPlacement.first, newPlacement.second);
            }
          }
        }

        // Still block the locations that are hit too early
        if (config.state.gen < params.startGen) {
          for (int s = 0; s < catalysts.size(); s++) {
            LifeState hitLocations = newcells.Convolve(catalystReactionMasks[s]);
            masks[s].Join(hitLocations);
          }
        }
      }

      if (config.count == params.numCatalysts) {
        bool allRecovered = true;
        for (int i = 0; i < config.count; i++) {
          if (!hasRecovered[i]) {
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
}
