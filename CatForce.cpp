// CatForce - Catalyst search utility based on LifeAPI using brute force.
// Written by Michael Simkin 2015
#include "LifeAPI.h"
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <stack>

const bool DEBUG_OUTPUT = false;
const int MAX_CATALYSTS = 5;
const int REQUIRED_LOOKAHEAD = 5;
const int LIGHTSPEED_LOOKAHEAD = 20;

void UpdateCounts(LifeState &__restrict__ state, LifeState &__restrict__ out1, LifeState &__restrict__ out2, LifeState &__restrict__ outMore) {
  LifeState bit3(false), bit2(false), bit1(false), bit0(false);
  state.CountNeighbourhood(bit3, bit2, bit1, bit0);
  out1 |= ~state & ~bit3 & ~bit2 & ~bit1 & bit0;
  out2 |= ~state & ~bit3 & ~bit2 & bit1 & ~bit0;
  outMore |= ~state & (bit3 | bit2 | (bit1 & bit0));
}

void SetCounts(LifeState &__restrict__ state, LifeState &__restrict__ out1, LifeState &__restrict__ out2, LifeState &__restrict__ outMore) {
  LifeState bit3(false), bit2(false), bit1(false), bit0(false);
  state.CountNeighbourhood(bit3, bit2, bit1, bit0);
  out1 = ~state & ~bit3 & ~bit2 & ~bit1 & bit0;
  out2 = ~state & ~bit3 & ~bit2 & bit1 & ~bit0;
  outMore = ~state & (bit3 | bit2 | (bit1 & bit0));
}


void split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
}

std::vector<std::string> splitwhitespace(std::string const &input) {
    std::istringstream buffer(input);
    std::vector<std::string> ret((std::istream_iterator<std::string>(buffer)),
                                  std::istream_iterator<std::string>());
    return ret;
}

class SearchParams {
public:
  unsigned maxGen;
  unsigned numCatalysts;
  unsigned numTransparent;
  unsigned numLimited;
  unsigned stableInterval;
  std::string pat;
  int xPat;
  int yPat;
  unsigned startGen;
  unsigned lastGen;
  std::string outputFile;
  std::string fullReportFile;
  int searchArea[4];
  int maxW;
  int maxH;
  StaticSymmetry symmetry;
  std::vector<SymmetryTransform> symmetryChain;

  int maxCatSize;

  std::string alsoRequired;
  std::pair<int, int> alsoRequiredXY;

  // Filtering parameters
  int stopAfterCatsDestroyed;
  int maxJunk;
  int matchSurvive;

  int offsetArea[4]{};

  bool useCollisionMasks;

  SearchParams() {
    maxGen = 250;
    numCatalysts = 2;
    numTransparent = 100;
    numLimited = 100;
    stableInterval = 15;
    pat = "";
    searchArea[0] = -30;
    searchArea[1] = -30;
    searchArea[2] = 60;
    searchArea[3] = 60;
    offsetArea[0] = -15;
    offsetArea[1] = -15;
    offsetArea[2] = 30;
    offsetArea[3] = 30;
    xPat = 0;
    yPat = 0;
    startGen = 0;
    lastGen = 100;
    outputFile = "results.rle";
    fullReportFile = "";
    maxW = -1;
    maxH = -1;
    symmetry = StaticSymmetry::C1;
    symmetryChain = {};
    maxCatSize = -1;
    alsoRequired = "";
    alsoRequiredXY = {0, 0};
    stopAfterCatsDestroyed = -1;
    maxJunk = -1;
    matchSurvive = -1;
    useCollisionMasks = true;
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
  std::string antirequiredRLE;
  std::pair<int, int> antirequiredXY;
  std::string locusRLE;
  std::pair<int, int> locusXY;
  std::string contactRLE;
  std::pair<int, int> contactXY;
  bool transparent;
  bool limited;
  bool mustInclude;
  bool checkRecovery;
  bool checkReaction;
  bool canSmother;
  bool canRock;
  bool sacrificial;
  bool periodic;
  bool fixed;
  unsigned fixedGen;

  explicit CatalystInput(std::string &line) {
    std::vector<std::string> elems = splitwhitespace(line);

    if (elems.size() < 6) {
      std::cout << "The line " << line << "is invalid" << std::endl;
      std::cout << "Format: cat <rle> <absense interval> <centerX> <centerY> "
                   "<symm Type | + / x *>"
                << std::endl;
      exit(1);
    }

    rle = elems[1];
    maxDisappear = atoi(elems[2].c_str());
    if (maxDisappear > 100) {
      std::cout << line << std::endl;
      exit(1);
    }
    centerX = atoi(elems[3].c_str());
    centerY = atoi(elems[4].c_str());
    symmType = elems[5].at(0);

    transparent = false;
    limited = false;
    mustInclude = false;
    checkRecovery = false;
    checkReaction = false;
    canSmother = false;
    canRock = false;
    sacrificial = false;
    periodic = false;
    fixed = false;
    fixedGen = 0;

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
      } else if (elems[argi] == "antirequired") {
        antirequiredRLE = elems[argi + 1];
        antirequiredXY = std::make_pair(atoi(elems[argi + 2].c_str()), atoi(elems[argi + 3].c_str()));
        argi += 4;
      } else if (elems[argi] == "locus") {
        locusRLE = elems[argi + 1];
        locusXY = std::make_pair(atoi(elems[argi + 2].c_str()), atoi(elems[argi + 3].c_str()));
        argi += 4;
      } else if (elems[argi] == "contact") {
        contactRLE = elems[argi + 1];
        contactXY = std::make_pair(atoi(elems[argi + 2].c_str()), atoi(elems[argi + 3].c_str()));
        argi += 4;
      } else if (elems[argi] == "transparent") {
        transparent = true;
        argi += 1;
      } else if (elems[argi] == "limited") {
        limited = true;
        argi += 1;
      } else if (elems[argi] == "mustinclude" || elems[argi] == "must-include") {
        mustInclude = true;
        argi += 1;
      } else if (elems[argi] == "check-recovery") {
        checkRecovery = true;
        argi += 1;
      } else if (elems[argi] == "check-reaction") {
        checkReaction = true;
        argi += 1;
      } else if (elems[argi] == "can-smother") {
        canSmother = true;
        argi += 1;
      } else if (elems[argi] == "can-rock") {
        canRock = true;
        argi += 1;
      } else if (elems[argi] == "sacrificial") {
        sacrificial = true;
        argi += 1;
      } else if (elems[argi] == "period") {
        periodic = true;
        argi += 2;
      } else if (elems[argi] == "fixed") {
        fixed = true;
        argi += 1;
      } else if (elems[argi] == "fixed-gen") {
        fixed = true;
        fixedGen = atoi(elems[argi + 1].c_str());
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

enum FilterType {
  ANDFILTER,
  ORFILTER,
  MATCHFILTER,
};

class FilterInput {
public:
  std::string rle;
  int x;
  int y;
  int gen;
  std::pair<int, int> range;
  FilterType type;
  char sym;

  explicit FilterInput(std::string &line) {
    std::vector<std::string> elems = splitwhitespace(line);

    std::vector<std::string> rangeElems;
    split(elems[1], '-', rangeElems);

    if (rangeElems.size() == 1) {
      gen = atoi(elems[1].c_str());
      range = {-1, -1};
    } else {
      unsigned minGen = atoi(rangeElems[0].c_str());
      unsigned maxGen = atoi(rangeElems[1].c_str());

      gen = -1;
      range = {minGen, maxGen};
    }

    rle = elems[2];

    std::string filter = "filter";
    std::string andfilter = "andfilter";
    std::string orfilter = "orfilter";
    std::string matchfilter = "match";

    if (elems[0] == andfilter) {
      x = atoi(elems[3].c_str());
      y = atoi(elems[4].c_str());
      type = ANDFILTER;
    } else if (elems[0] == orfilter) {
      x = atoi(elems[3].c_str());
      y = atoi(elems[4].c_str());
      type = ORFILTER;
    } else if (elems[0] == matchfilter) {
      x = 0;
      y = 0;
      type = MATCHFILTER;
      if (elems.size() > 3)
        sym = elems[3].at(0);
      else
        sym = '*';
    } else {
      x = atoi(elems[3].c_str());
      y = atoi(elems[4].c_str());
      type = ANDFILTER;
    }
  }
};

inline std::vector<SymmetryTransform> AllTransforms() {
  return {Identity,
          ReflectAcrossXEven,
          ReflectAcrossX,
          ReflectAcrossYEven,
          ReflectAcrossY,
          Rotate90Even,
          Rotate90,
          Rotate270Even,
          Rotate270,
          Rotate180OddBoth,
          Rotate180EvenHorizontal,
          Rotate180EvenVertical,
          Rotate180EvenBoth,
          ReflectAcrossYeqX,
          ReflectAcrossYeqNegX,
          ReflectAcrossYeqNegXP1};
}

inline std::vector<SymmetryTransform> SymmetryGroupFromEnum(const StaticSymmetry sym) {
  switch (sym) {
  case StaticSymmetry::C1:
    return {Identity};
  case StaticSymmetry::D2AcrossX:
    return {Identity, ReflectAcrossX};
    // vertical/horizontal here refer to box dimensions, NOT axis of reflection
  case StaticSymmetry::D2AcrossXEven:
    return {Identity, ReflectAcrossXEven};
  case StaticSymmetry::D2AcrossY:
    return {Identity, ReflectAcrossY};
  case StaticSymmetry::D2AcrossYEven:
    return {Identity, ReflectAcrossYEven};
  case StaticSymmetry::D2diagodd:
    return {Identity, ReflectAcrossYeqX};
  case StaticSymmetry::D2negdiagodd:
    return {Identity, ReflectAcrossYeqNegXP1};
  case StaticSymmetry::C2:
    return {Identity, Rotate180OddBoth};
  case StaticSymmetry::C2even:
    return {Identity, Rotate180EvenBoth};
  case StaticSymmetry::C2horizontaleven:
    return {Identity, Rotate180EvenHorizontal};
  case StaticSymmetry::C2verticaleven:
    return {Identity, Rotate180EvenVertical};
  case StaticSymmetry::C4:
    return {Identity, Rotate90, Rotate180OddBoth, Rotate270};
  case StaticSymmetry::C4even:
    return {Identity, Rotate90Even, Rotate180EvenBoth, Rotate270Even};
  case StaticSymmetry::D4:
    return {Identity, ReflectAcrossX, Rotate180OddBoth, ReflectAcrossY};
  case StaticSymmetry::D4even:
    return {Identity, ReflectAcrossXEven, Rotate180EvenBoth,
            ReflectAcrossYEven};
  case StaticSymmetry::D4horizontaleven:
    return {Identity, ReflectAcrossYEven, Rotate180EvenHorizontal,
            ReflectAcrossX};
  case StaticSymmetry::D4verticaleven:
    return {Identity, ReflectAcrossXEven, Rotate180EvenVertical,
            ReflectAcrossY};
  case StaticSymmetry::D4diag:
    return {Identity, ReflectAcrossYeqX, Rotate180OddBoth,
            ReflectAcrossYeqNegXP1};
  case StaticSymmetry::D4diageven:
    return {Identity, ReflectAcrossYeqX, Rotate180EvenBoth,
            ReflectAcrossYeqNegX};
  case StaticSymmetry::D8:
    return {Identity,       ReflectAcrossX,         ReflectAcrossYeqX,
            ReflectAcrossY, ReflectAcrossYeqNegXP1, Rotate90,
            Rotate270,      Rotate180OddBoth};
  case StaticSymmetry::D8even:
    return {Identity,           ReflectAcrossXEven,   ReflectAcrossYeqX,
            ReflectAcrossYEven, ReflectAcrossYeqNegX, Rotate90Even,
            Rotate270Even,      Rotate180EvenBoth};
  }
}

inline std::vector<SymmetryTransform> SymmetryChainFromEnum(const StaticSymmetry sym) {
  switch (sym) {
  case StaticSymmetry::C1:
    return {};
  case StaticSymmetry::D2AcrossY:
    return {ReflectAcrossY};
  case StaticSymmetry::D2AcrossYEven:
    return {ReflectAcrossYEven};
  case StaticSymmetry::D2AcrossX:
    return {ReflectAcrossX};
  case StaticSymmetry::D2AcrossXEven:
    return {ReflectAcrossXEven};
  case StaticSymmetry::D2diagodd:
    return {ReflectAcrossYeqX};
  case StaticSymmetry::D2negdiagodd:
    return {ReflectAcrossYeqNegXP1};
  case StaticSymmetry::C2:
    return {Rotate180OddBoth};
  case StaticSymmetry::C2even:
    return {Rotate180EvenBoth};
  case StaticSymmetry::C2horizontaleven:
    return {Rotate180EvenHorizontal};
  case StaticSymmetry::C2verticaleven:
    return {Rotate180EvenVertical};
  case StaticSymmetry::C4:
    return {Rotate90, Rotate180OddBoth};
  case StaticSymmetry::C4even:
    return {Rotate90Even, Rotate180EvenBoth};
  case StaticSymmetry::D4:
    return {ReflectAcrossX, ReflectAcrossY};
  case StaticSymmetry::D4even:
    return {ReflectAcrossXEven, ReflectAcrossYEven};
  case StaticSymmetry::D4horizontaleven:
    return {ReflectAcrossYEven, ReflectAcrossX};
  case StaticSymmetry::D4verticaleven:
    return {ReflectAcrossXEven, ReflectAcrossY};
  case StaticSymmetry::D4diag:
    return {ReflectAcrossYeqX, ReflectAcrossYeqNegXP1};
  case StaticSymmetry::D4diageven:
    return {ReflectAcrossYeqX, ReflectAcrossYeqNegX};
  case StaticSymmetry::D8:
    return {Rotate90, Rotate180OddBoth, ReflectAcrossYeqX};
  case StaticSymmetry::D8even:
    return {Rotate90Even, Rotate180EvenBoth, ReflectAcrossYeqX};
  }
}

const LifeState C1Domain = LifeState::Parse(
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$"
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$"
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$"
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o!");

const LifeState D2YDomain = LifeState::Parse(
        "32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$"
        "32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$"
        "32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$"
        "32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o!");

const LifeState D2XDomain = LifeState::Parse(
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$"
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o!");

const LifeState D2DiagDomain = LifeState::Parse(
    "o$2o$3o$4o$5o$6o$7o$8o$9o$10o$11o$12o$13o$14o$15o$16o$17o$18o$19o$20o$"
    "21o$22o$23o$24o$25o$26o$27o$28o$29o$30o$31o$32o$33o$34o$35o$36o$37o$"
    "38o$39o$40o$41o$42o$43o$44o$45o$46o$47o$48o$49o$50o$51o$52o$53o$54o$"
    "55o$56o$57o$58o$59o$60o$61o$62o$63o$64o!");

const LifeState D2NegDiagDomain = LifeState::Parse(
    "64o$63o$62o$61o$60o$59o$58o$57o$56o$55o$54o$53o$52o$51o$50o$49o$48o$"
    "47o$46o$45o$44o$43o$42o$41o$40o$39o$38o$37o$36o$35o$34o$33o$32o$31o$"
    "30o$29o$28o$27o$26o$25o$24o$23o$22o$21o$20o$19o$18o$17o$16o$15o$14o$"
    "13o$12o$11o$10o$9o$8o$7o$6o$5o$4o$3o$2o$o!");

const LifeState C2Domain = LifeState::Parse(
    "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$"
    "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o!");

const LifeState C4Domain = LifeState::Parse(
    "32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$"
    "32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o!");

const LifeState D4DiagDomain = LifeState::Parse(
    "o$2o$3o$4o$5o$6o$7o$8o$9o$10o$11o$12o$13o$14o$15o$16o$17o$18o$19o$20o$"
    "21o$22o$23o$24o$25o$26o$27o$28o$29o$30o$31o$32o$32o$31o$30o$29o$28o$"
    "27o$26o$25o$24o$23o$22o$21o$20o$19o$18o$17o$16o$15o$14o$13o$12o$11o$"
    "10o$9o$8o$7o$6o$5o$4o$3o$2o$o!");

const LifeState D8Domain = LifeState::Parse(
        "o$2o$3o$4o$5o$6o$7o$8o$9o$10o$11o$12o$13o$14o$15o$16o$17o$18o$19o$20o$"
        "21o$22o$23o$24o$25o$26o$27o$28o$29o$30o$31o$32o!");

LifeState FundamentalDomainFast(const StaticSymmetry sym) {
  switch (sym) {
  case StaticSymmetry::C1:
    return C1Domain;
  case StaticSymmetry::D2AcrossY:
  case StaticSymmetry::D2AcrossYEven:
    return D2YDomain;
  case StaticSymmetry::D2AcrossX:
  case StaticSymmetry::D2AcrossXEven:
    return D2XDomain;
  case StaticSymmetry::D2diagodd:
    return D2DiagDomain;
  case StaticSymmetry::D2negdiagodd:
    return D2NegDiagDomain;
  case StaticSymmetry::C2:
  case StaticSymmetry::C2even:
  case StaticSymmetry::C2horizontaleven:
  case StaticSymmetry::C2verticaleven:
    return C2Domain;
  case StaticSymmetry::C4:
  case StaticSymmetry::C4even:
  case StaticSymmetry::D4:
  case StaticSymmetry::D4even:
  case StaticSymmetry::D4horizontaleven:
  case StaticSymmetry::D4verticaleven:
    return C4Domain;
  case StaticSymmetry::D4diag:
  case StaticSymmetry::D4diageven:
    return D4DiagDomain;
  case StaticSymmetry::D8:
  case StaticSymmetry::D8even:
    return D8Domain;
  }
}

LifeState FundamentalDomain(const StaticSymmetry sym) {
  switch (sym) {
  case StaticSymmetry::C1:
    return LifeState::Parse(
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$"
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$"
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$"
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o!");
  case StaticSymmetry::D2AcrossY:
  case StaticSymmetry::D2AcrossYEven:
    return LifeState::Parse(
        "32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$"
        "32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$"
        "32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$"
        "32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o!");
  case StaticSymmetry::D2AcrossX:
  case StaticSymmetry::D2AcrossXEven:
    return LifeState::Parse(
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$"
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o!");
  case StaticSymmetry::D2diagodd:
    return LifeState::Parse(
        "o$2o$3o$4o$5o$6o$7o$8o$9o$10o$11o$12o$13o$14o$15o$16o$17o$18o$19o$20o$"
        "21o$22o$23o$24o$25o$26o$27o$28o$29o$30o$31o$32o$33o$34o$35o$36o$37o$"
        "38o$39o$40o$41o$42o$43o$44o$45o$46o$47o$48o$49o$50o$51o$52o$53o$54o$"
        "55o$56o$57o$58o$59o$60o$61o$62o$63o$64o!");
  case StaticSymmetry::D2negdiagodd:
    return LifeState::Parse(
        "64o$63o$62o$61o$60o$59o$58o$57o$56o$55o$54o$53o$52o$51o$50o$49o$48o$"
        "47o$46o$45o$44o$43o$42o$41o$40o$39o$38o$37o$36o$35o$34o$33o$32o$31o$"
        "30o$29o$28o$27o$26o$25o$24o$23o$22o$21o$20o$19o$18o$17o$16o$15o$14o$"
        "13o$12o$11o$10o$9o$8o$7o$6o$5o$4o$3o$2o$o!");
  case StaticSymmetry::C2:
  case StaticSymmetry::C2even:
  case StaticSymmetry::C2horizontaleven:
  case StaticSymmetry::C2verticaleven:
    return LifeState::Parse(
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$"
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o!");
  case StaticSymmetry::C4:
  case StaticSymmetry::C4even:
  case StaticSymmetry::D4:
  case StaticSymmetry::D4even:
  case StaticSymmetry::D4horizontaleven:
  case StaticSymmetry::D4verticaleven:
    return LifeState::Parse(
        "32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$"
        "32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o!");
  case StaticSymmetry::D4diag:
  case StaticSymmetry::D4diageven:
    return LifeState::Parse(
        "o$2o$3o$4o$5o$6o$7o$8o$9o$10o$11o$12o$13o$14o$15o$16o$17o$18o$19o$20o$"
        "21o$22o$23o$24o$25o$26o$27o$28o$29o$30o$31o$32o$32o$31o$30o$29o$28o$"
        "27o$26o$25o$24o$23o$22o$21o$20o$19o$18o$17o$16o$15o$14o$13o$12o$11o$"
        "10o$9o$8o$7o$6o$5o$4o$3o$2o$o!");
  case StaticSymmetry::D8:
  case StaticSymmetry::D8even:
    return LifeState::Parse(
        "o$2o$3o$4o$5o$6o$7o$8o$9o$10o$11o$12o$13o$14o$15o$16o$17o$18o$19o$20o$"
        "21o$22o$23o$24o$25o$26o$27o$28o$29o$30o$31o$32o!");
  }
}



inline std::pair<int, int> CommuteTranslation(const SymmetryTransform sym,
                                              std::pair<int, int> vec) {
  int x = vec.first;
  int y = vec.second;
  switch (sym) {
  case Identity:
    return std::make_pair(x, y);
  case ReflectAcrossXEven:
    return std::make_pair(x, -y);
  case ReflectAcrossX:
    return std::make_pair(x, -y);
  case ReflectAcrossYEven:
    return std::make_pair(-x, y);
  case ReflectAcrossY:
    return std::make_pair(-x, y);
  case Rotate90Even:
    return std::make_pair(-y, x);
  case Rotate90:
    return std::make_pair(-y, x);
  case Rotate270Even:
    return std::make_pair(y, -x);
  case Rotate270:
    return std::make_pair(y, -x);
  case Rotate180OddBoth:
    return std::make_pair(-x, -y);
  case Rotate180EvenHorizontal:
    return std::make_pair(-x, -y);
  case Rotate180EvenVertical:
    return std::make_pair(-x, -y);
  case Rotate180EvenBoth:
    return std::make_pair(-x, -y);
  case ReflectAcrossYeqX:
    return std::make_pair(y, x);
  case ReflectAcrossYeqNegX:
    return std::make_pair(-y, -x);
  case ReflectAcrossYeqNegXP1:
    return std::make_pair(-y, -x);
  }
}

inline std::pair<int, int> HalveOffset(const StaticSymmetry sym,
                                       std::pair<int, int> vec) {
  switch (sym) {
  case C4: { // This is the center of rotation for offset rotation by 90
    int x = vec.first;
    int y = vec.second;
    int x2 = (x - y) / 2;
    int y2 = (x + y) / 2;
    int x3 = ((x2 + 16 + 32) % 32 - 16 + 64) % 64;
    int y3 = ((y2 + 16 + 32) % 32 - 16 + 64) % 64;
    return std::make_pair(x3, y3);
  }
  default: {
    int x = (((vec.first + 32) % 64 - 32) / 2 + 64) % 64;
    int y = (((vec.second + 32) % 64 - 32) / 2 + 64) % 64;
    return std::make_pair(x, y);
  }
  }
}

StaticSymmetry SymmetryFromString(const std::string &name) {
  std::string start = name.substr(0, 2);
  std::string rest = name.substr(2);
  if (start == "D2") {
    if (rest == "-" or rest == "vertical") {
      return StaticSymmetry::D2AcrossX;
    } else if (rest == "-even" or rest == "verticaleven") {
      return StaticSymmetry::D2AcrossXEven;
    } else if (rest == "|" or rest == "horizontal") {
      return StaticSymmetry::D2AcrossY;
    } else if (rest == "|even" or rest == "horizontaleven") {
      return StaticSymmetry::D2AcrossYEven;
    } else if (rest == "/" or rest == "/odd") {
      return StaticSymmetry::D2negdiagodd;
    } else if (rest == "\\" or rest == "\\odd") {
      return StaticSymmetry::D2diagodd;
    }
  } else if (start == "C2") {
    if (rest == "" or rest == "_1") {
      return StaticSymmetry::C2;
    } else if (rest == "even" or rest == "_4") {
      return StaticSymmetry::C2even;
    } else if (rest == "horizontaleven" or rest == "|even") {
      return StaticSymmetry::C2horizontaleven;
    } else if (rest == "verticaleven" or rest == "-even" or rest == "_2") {
      return StaticSymmetry::C2verticaleven;
    }
  } else if (start == "C4") {
    if (rest == "" or rest == "_1") {
      return StaticSymmetry::C4;
    } else if (rest == "even" or rest == "_4") {
      return StaticSymmetry::C4even;
    }
  } else if (start == "D4") {
    std::string evenOddInfo = rest.substr(1);
    if (rest[0] == '+' or (rest.size() > 1 and rest[1] == '+')) {
      if (evenOddInfo == "" or rest == "_+1") {
        return StaticSymmetry::D4;
      } else if (evenOddInfo == "even" or rest == "_+4") {
        return StaticSymmetry::D4even;
      } else if (evenOddInfo == "verticaleven" or evenOddInfo == "-even" or
                 rest == "_+2") {
        return StaticSymmetry::D4verticaleven;
      } else if (evenOddInfo == "horizontaleven" or evenOddInfo == "|even") {
        return StaticSymmetry::D4horizontaleven;
      }
    } else if (rest[0] == 'x' or (rest.size() > 1 and rest[1] == 'x')) {
      if (evenOddInfo == "" or rest == "_x1") {
        return StaticSymmetry::D4diag;
      } else if (evenOddInfo == "even" or rest == "_x4") {
        return StaticSymmetry::D4diageven;
      }
    }
  } else if (start == "D8") {
    if (rest == "" or rest == "_1") {
      return StaticSymmetry::D8;
    } else if (rest == "even" or rest == "_4") {
      return StaticSymmetry::D8even;
    }
  }
  return StaticSymmetry::C1;
}

std::string SymmetryToString(StaticSymmetry sym) {
  switch (sym) {
  case StaticSymmetry::C1:
    return "C1";
  case StaticSymmetry::D2AcrossX:
    return "D2-";
  case StaticSymmetry::D2AcrossXEven:
    return "D2-even";
  case StaticSymmetry::D2AcrossY:
    return "D2|";
  case StaticSymmetry::D2AcrossYEven:
    return "D2|even";
  case StaticSymmetry::D2diagodd:
    return "D2\\";
  case StaticSymmetry::D2negdiagodd:
    return "D2/";
  case StaticSymmetry::C2:
    return "C2";
  case StaticSymmetry::C2even:
    return "C2even";
  case StaticSymmetry::C2horizontaleven:
    return "C2|even";
  case StaticSymmetry::C2verticaleven:
    return "C2-even";
  case StaticSymmetry::C4:
    return "C4";
  case StaticSymmetry::C4even:
    return "C4even";
  case StaticSymmetry::D4:
    return "D4+";
  case StaticSymmetry::D4even:
    return "D4+even";
  case StaticSymmetry::D4horizontaleven:
    return "D4+|even";
  case StaticSymmetry::D4verticaleven:
    return "D4+-even";
  case StaticSymmetry::D4diag:
    return "D4x";
  case StaticSymmetry::D4diageven:
    return "D4xeven";
  case StaticSymmetry::D8:
    return "D8";
  case StaticSymmetry::D8even:
    return "D8even";
  }
}

std::vector<SymmetryTransform> CharToTransforms(char ch) {
  switch (ch) {
  case '.':
    return SymmetryGroupFromEnum(StaticSymmetry::C1);
  case '|':
    return SymmetryGroupFromEnum(StaticSymmetry::D2AcrossY);
  case '-':
    return SymmetryGroupFromEnum(StaticSymmetry::D2AcrossX);
  case '\\':
    return SymmetryGroupFromEnum(StaticSymmetry::D2diagodd);
  case '/':
    return SymmetryGroupFromEnum(StaticSymmetry::D2negdiagodd);
  case '+':
  case '@':
    return SymmetryGroupFromEnum(StaticSymmetry::C4);
  case 'x':
    return {Identity, Rotate90, ReflectAcrossX, ReflectAcrossYeqX};
  case '*':
    return SymmetryGroupFromEnum(StaticSymmetry::D8);
  default:
    return SymmetryGroupFromEnum(StaticSymmetry::C1);
  }
}

void ReadParams(const std::string &fname, std::vector<CatalystInput> &catalysts,
                std::vector<FilterInput> &filters,
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
  std::string numLimited = "num-limited";
  std::string stable = "stable-interval";
  std::string area = "search-area";
  std::string pat = "pat";
  std::string outputFile = "output";
  std::string filter = "filter";
  std::string andfilter = "andfilter";
  std::string orfilter = "orfilter";
  std::string matchfilter = "match";
  std::string maxWH = "fit-in-width-height";
  std::string maxCatSize = "max-category-size";
  std::string fullReport = "full-report";

  std::string symmetry = "symmetry";
  std::string alsoRequired = "also-required";
  std::string stopAfterCatsDestroyed = "stop-after-cats-destroyed";
  std::string maxJunk = "max-junk";
  std::string matchSurvive = "match-survive";

  std::string offsetArea = "offset-area";
  std::string useCollisionMasks = "use-collision-masks";

  std::string line;
  while (std::getline(infile, line)) {
    std::vector<std::string> elems = splitwhitespace(line);

    if (elems.size() < 2)
      continue;

    if (elems[0] == Cat) {
      catalysts.emplace_back(line);
    } else if (elems[0] == maxGen) {
      params.maxGen = atoi(elems[1].c_str());
    } else if (elems[0] == numCat) {
      params.numCatalysts = atoi(elems[1].c_str());
      if(params.numCatalysts > MAX_CATALYSTS) {
        std::cout << "num-catalyst larger than MAX_CATALYSTS" << std::endl;
        exit(1);
      }
    } else if (elems[0] == numTransp) {
      params.numTransparent = atoi(elems[1].c_str());
    } else if (elems[0] == numLimited) {
      params.numLimited = atoi(elems[1].c_str());
    } else if (elems[0] == stable) {
      params.stableInterval = atoi(elems[1].c_str());
    } else if (elems[0] == pat) {
      params.pat = elems[1];

      if (elems.size() > 3) {
        params.xPat = atoi(elems[2].c_str());
        params.yPat = atoi(elems[3].c_str());
      }
    } else if (elems[0] == area) {
      params.searchArea[0] = atoi(elems[1].c_str());
      params.searchArea[1] = atoi(elems[2].c_str());
      params.searchArea[2] = atoi(elems[3].c_str());
      params.searchArea[3] = atoi(elems[4].c_str());
    } else if (elems[0] == startGen) {
      params.startGen = atoi(elems[1].c_str());
    } else if (elems[0] == lastGen) {
      params.lastGen = atoi(elems[1].c_str());
    } else if (elems[0] == outputFile) {
      params.outputFile = elems[1];

      for (unsigned i = 2; i < elems.size(); i++) {
        params.outputFile.append(" ");
        params.outputFile.append(elems[i]);
      }
    } else if (elems[0] == fullReport) {
      params.fullReportFile = elems[1];

      for (unsigned i = 2; i < elems.size(); i++) {
        params.fullReportFile.append(" ");
        params.fullReportFile.append(elems[i]);
      }
    } else if (elems[0] == filter || elems[0] == orfilter || elems[0] == andfilter || elems[0] == matchfilter) {
      filters.emplace_back(line);
    } else if (elems[0] == maxWH) {
      params.maxW = atoi(elems[1].c_str());
      params.maxH = atoi(elems[2].c_str());
    } else if (elems[0] == maxCatSize) {
      params.maxCatSize = atoi(elems[1].c_str());
    } else if (elems[0] == symmetry) {
      // std::string symmetryString = "";
      // // reverse-compatibility reasons.
      // if (elems[1] == "horizontal") {
      //   symmetryString = "D2|odd";
      // } else if (elems[1] == "horizontaleven") {
      //   symmetryString = "D2|even";
      // } else if (elems[1] == "diagonal") {
      //   symmetryString =
      //       "D2/"; // I think this was the way that it worked before?
      // } else if (elems[1] == "rotate180") {
      //   symmetryString = "C2";
      // } else if (elems[1] == "rotate180evenx") {
      //   symmetryString = "C2horizontaleven";
      // } else if (elems[1] == "rotate180evenboth") {
      //   symmetryString = "C2evenboth";
      // } else {
      //   symmetryString = elems[1];
      // }

      // params.symmetry = SymmetryFromString(symmetryString);
      // params.symmetryChain = SymmetryChainFromEnum(params.symmetry);

      std::cout << "symmetry option forbidden in offsets branch" << std::endl;
      exit(1);
    } else if (elems[0] == alsoRequired) {
      params.alsoRequired = elems[1].c_str();
      params.alsoRequiredXY = std::make_pair(atoi(elems[2].c_str()), atoi(elems[3].c_str()));
    } else if (elems[0] == stopAfterCatsDestroyed){
      params.stopAfterCatsDestroyed = atoi(elems[1].c_str());
    } else if (elems[0] == maxJunk){
      params.maxJunk = atoi(elems[1].c_str());
    } else if (elems[0] == offsetArea) {
      params.offsetArea[0] = atoi(elems[1].c_str());
      params.offsetArea[1] = atoi(elems[2].c_str());
      params.offsetArea[2] = atoi(elems[3].c_str());
      params.offsetArea[3] = atoi(elems[4].c_str());
    } else if (elems[0] == matchSurvive){
      params.matchSurvive = atoi(elems[1].c_str());
    } else if (elems[0] == useCollisionMasks){
      params.useCollisionMasks = atoi(elems[1].c_str());
    } else {
      if(std::isalpha(elems[0][0])) {
        std::cout << "Unknown input parameter: " << elems[0] << std::endl;
        exit(1);
      }
    }
  }
}

class CatalystData {
public:
  LifeState state;

  LifeState state1;
  LifeState state2;
  LifeState stateMore;
  LifeState statezoi;

  LifeState reaction1;
  LifeState reaction2;
  LifeState reactionMore;

  // bool hasLocusReactionPop;
  bool hasLocusReactionPop1;
  bool hasLocusReactionPop2;
  // bool hasLocusReactionPopMore;

  // unsigned locusReactionPop;
  unsigned locusReactionPop1;
  unsigned locusReactionPop2;
  // unsigned locusReactionPopMore;

  LifeState locusReactionMask;
  LifeState locusReactionMask1;  // Positions that react with a cell with 1 neighbour at the origin etc.
  LifeState locusReactionMask2;
  LifeState locusReactionMaskMore;

  LifeState locusAvoidMask;
  LifeState locusAvoidMask1;
  LifeState locusAvoidMask2;
  LifeState locusAvoidMaskMore;

  // p2 HACK
  LifeState offReactionMask;
  LifeState off;
  LifeState off1;
  LifeState off2;
  LifeState offMore;
  LifeState offReactionMask1;
  LifeState offReactionMask2;
  LifeState offReactionMaskMore;

  LifeState required;

  LifeTarget target;

  bool hasRequired;
  bool hasLocus;

  std::vector<LifeTarget> forbidden;

  unsigned maxDisappear;
  bool transparent;
  bool limited;
  bool mustInclude;
  bool checkRecovery;
  bool checkReaction;
  bool canSmother;
  bool canRock;
  bool sacrificial;
  bool periodic;
  unsigned phase;
  bool fixed;
  unsigned fixedGen;

  static std::vector<CatalystData> FromInput(CatalystInput &input);
};

std::vector<CatalystData> CatalystData::FromInput(CatalystInput &input) {
  std::vector<SymmetryTransform> trans = CharToTransforms(input.symmType);

  const char *rle = input.rle.c_str();

  std::vector<CatalystData> results;

  for (int phase = 0; phase < 2; phase++) {
  for (auto &tran : trans) {
    LifeState pat = LifeState::Parse(rle, input.centerX, input.centerY, tran);
    LifeState phasepat = pat;
    phasepat.Step(phase);

    CatalystData result;
    result.state = phasepat;
    result.statezoi = phasepat.ZOI();

    LifeState gen1 = pat;
    gen1.Step();
    result.target = LifeTarget();
    result.target.wanted = pat & gen1;
    result.target.unwanted = (pat | gen1).ZOI() & ~(pat | gen1);

    result.hasLocus = false;

    result.required = LifeState();

    if (input.requiredRLE != "") {
      result.hasRequired = true;
      result.required |= LifeState::Parse(input.requiredRLE.c_str(),
                                          input.requiredXY.first,
                                          input.requiredXY.second, tran);
    }

    if (input.antirequiredRLE != "") {
      result.hasRequired = true;
      result.required |= LifeState::Parse(input.antirequiredRLE.c_str(),
                                          input.antirequiredXY.first,
                                          input.antirequiredXY.second, tran);
    }

    result.required &= ~(pat ^ gen1);

    LifeState locus;
    if (input.locusRLE != "") {
      result.hasLocus = true;
      locus = LifeState::Parse(input.locusRLE.c_str(),
                                      input.locusXY.first,
                                      input.locusXY.second, tran);
    } else {
      locus = pat;
    }

    LifeState shell = pat.ZOI().Shell();

    LifeState contact;
    if (input.contactRLE != "") {
      contact = LifeState::Parse(input.contactRLE.c_str(),
                                 input.contactXY.first,
                                 input.contactXY.second, tran);
    } else {
      LifeState locusZOI = locus.ZOI();
      LifeState nonLocusZOI = (pat & ~locus).ZOI();
      contact = locusZOI & ~nonLocusZOI & shell & ~result.required;
    }

    {
      result.locusAvoidMask = shell & ~contact;
      result.locusAvoidMask.Transform(Rotate180OddBoth);
      result.locusAvoidMask.RecalculateMinMax();

      result.locusReactionMask = shell & contact;
      result.locusReactionMask.Transform(Rotate180OddBoth);
      result.locusReactionMask.RecalculateMinMax();
      result.locusReactionMask &= ~result.locusAvoidMask;

      // result.locusAvoidMask = result.locusAvoidMask.Shell();

      LifeState state1, state2, stateMore;
      UpdateCounts(pat, state1, state2, stateMore);

      result.state1 = state1;
      result.state2 = state2;
      result.stateMore = stateMore;

      result.reaction1 = state1;
      result.reaction1.Transform(Rotate180OddBoth);
      result.reaction2 = state2;
      result.reaction2.Transform(Rotate180OddBoth);
      result.reactionMore = stateMore;
      result.reactionMore.Transform(Rotate180OddBoth);

      result.locusReactionMask1 = state1 & contact;
      result.locusReactionMask1.Transform(Rotate180OddBoth);
      result.locusReactionMask1.RecalculateMinMax();

      result.locusReactionMask2 = state2 & contact;
      result.locusReactionMask2.Transform(Rotate180OddBoth);
      result.locusReactionMask2.RecalculateMinMax();

      // result.locusReactionMaskMore = stateMore & shell & ~nonLocusZOI;
      // result.locusReactionMaskMore.Transform(Rotate180OddBoth);
      // result.locusReactionMaskMore.RecalculateMinMax();

      result.locusAvoidMask1 = state1 & shell & ~contact;
      result.locusAvoidMask1.Transform(Rotate180OddBoth);
      result.locusAvoidMask1.RecalculateMinMax();

      result.locusAvoidMask2 = state2 & shell & ~contact;
      result.locusAvoidMask2.Transform(Rotate180OddBoth);
      result.locusAvoidMask2.RecalculateMinMax();

      // result.locusAvoidMaskMore = stateMore & shell & nonLocusZOI;
      // result.locusAvoidMaskMore.Transform(Rotate180OddBoth);
      // result.locusAvoidMaskMore.RecalculateMinMax();

      // result.locusReactionPop = result.locusReactionMask.GetPop();
      result.locusReactionPop1 = result.locusReactionMask1.GetPop();
      result.locusReactionPop2 = result.locusReactionMask2.GetPop();
      // result.locusReactionPopMore = result.locusReactionMaskMore.GetPop();

      // result.hasLocusReactionPop = result.locusReactionPop > 0;
      result.hasLocusReactionPop1 = result.locusReactionPop1 > 0;
      result.hasLocusReactionPop2 = result.locusReactionPop2 > 0;
      // result.hasLocusReactionPopMore = result.locusReactionPopMore > 0;

      // p2 HACK
      result.offReactionMask = gen1.ZOI().Shell();
      result.offReactionMask.Transform(Rotate180OddBoth);
      result.offReactionMask.RecalculateMinMax();

      LifeState off1, off2, offMore;
      UpdateCounts(gen1, off1, off2, offMore);

      result.off = gen1;
      result.off1 = off1;
      result.off2 = off2;
      result.offMore = offMore;

      result.offReactionMask1 = off1;
      result.offReactionMask1.Transform(Rotate180OddBoth);
      result.offReactionMask2 = off2;
      result.offReactionMask2.Transform(Rotate180OddBoth);
      result.offReactionMaskMore = offMore;
      result.offMore.Transform(Rotate180OddBoth);
    }

    for (unsigned k = 0; k < input.forbiddenRLE.size(); k++) {
      result.forbidden.push_back(LifeTarget::Parse(input.forbiddenRLE[k].c_str(),
                                                   input.forbiddenXY[k].first,
                                                   input.forbiddenXY[k].second, tran));
    }

    result.maxDisappear = input.maxDisappear;
    result.transparent = input.transparent;
    result.limited = input.limited;
    result.mustInclude = input.mustInclude;
    result.checkRecovery = input.checkRecovery;
    result.checkReaction = input.checkReaction;
    result.canSmother = input.canSmother;
    result.canRock = input.canRock;
    result.sacrificial = input.sacrificial;
    result.periodic = input.periodic;
    result.phase = phase;
    result.fixed = input.fixed;
    result.fixedGen = input.fixedGen;

    if(input.fixed) {
      // We flip it back to avoid a convolve later
      result.locusReactionMask.Transform(Rotate180OddBoth);
      result.locusReactionMask.RecalculateMinMax();
    }

    results.push_back(result);
  }
  if(!input.periodic || input.fixed)
    break;
  }
  return results;
}

class FilterData {
public:
  LifeTarget target;

  std::vector<LifeTarget> transformedTargets;

  int gen;
  std::pair<int, int> range;
  FilterType type;

  static FilterData FromInput(FilterInput &input);
};

FilterData FilterData::FromInput(FilterInput &input) {
  FilterData result;

  result.target = LifeTarget::Parse(input.rle.c_str(), input.x, input.y);
  result.gen = input.gen;
  result.range = input.range;
  result.type = input.type;

  std::vector<SymmetryTransform> transforms = CharToTransforms(input.sym);
  for (auto trans : transforms) {
    LifeTarget transformed = result.target;
    transformed.Transform(trans);
    result.transformedTargets.push_back(transformed);
  }

  return result;
}


struct Configuration {
  LifeState startingCatalysts;
  unsigned count;
  unsigned transparentCount;
  unsigned limitedCount;
  unsigned mustIncludeCount;
  std::array<unsigned, MAX_CATALYSTS> curx;
  std::array<unsigned, MAX_CATALYSTS> cury;
  std::array<unsigned, MAX_CATALYSTS> curs;

  StaticSymmetry symmetry;
  std::pair<int, int> symmetryOffset;
};

// Fix a, what positions of b causes a collision?
LifeState CollisionMask(const CatalystData &a, const CatalystData &b) {
  LifeState bFlipped = b.state;
  bFlipped.Transform(Rotate180OddBoth);

  // Assumes both are still (so will break with periodic)
  LifeState possibleReactionMask =
      (a.state.ZOI().Convolve(bFlipped | b.reactionMore))
    | (a.state1.Convolve(b.reaction2))
    | (a.state2.Convolve(b.reaction1));

  if(a.periodic || b.periodic) {
    LifeState bOffFlipped = b.off;
    bOffFlipped.Transform(Rotate180OddBoth);

    possibleReactionMask |=
      (a.off.ZOI().Convolve(bOffFlipped | b.offReactionMaskMore))
      | (a.off1.Convolve(b.offReactionMask2))
      | (a.off2.Convolve(b.offReactionMask1));
  }

  return possibleReactionMask;
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
    categoryKey = catalystRemoved;
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

    LifeState result = afterCatalyst ^ catalysts;
    result.gen = firstGenSurvive;

    result.Step(maxgen - result.gen);
    uint64_t hash = result.GetHash();

    result = afterCatalyst ^ catalysts;
    result.gen = firstGenSurvive;

    for (auto & category: categories) {
      if (category->BelongsTo(result, hash)) {
          SearchResult r(init, conf, firstGenSurvive, genSurvive);
          category->Add(r);
          return;
      }
    }

    LifeState categoryKey = afterCatalyst ^ catalysts;
    categoryKey.gen = firstGenSurvive;

    SearchResult r(init, conf, firstGenSurvive, genSurvive);
    categories.push_back(new Category(categoryKey, r, catDelta, maxgen));
  }

  std::string CategoriesRLE(int maxCatSize) {
    std::stringstream ss;
    for (auto & category: categories) {
      ss << category->RLE(maxCatSize);
    }

    return ss.str();
  }
};

inline std::pair<int, int> PerpComponent(SymmetryTransform transf,
                                         std::pair<int, int> offset) {
  switch (transf) {
  case ReflectAcrossX:
    return std::make_pair(0, offset.second);
  case ReflectAcrossY:
    return std::make_pair(offset.first, 0);
  case ReflectAcrossYeqX: {
    int x = (offset.first + 32) % 64 - 32;
    int y = (offset.second + 32) % 64 - 32;
    return std::make_pair(((x - y + 128) / 2) % 64,
                          ((-x + y + 128) / 2) % 64);
  }
  case ReflectAcrossYeqNegXP1: {
    int x = (offset.first + 32) % 64 - 32;
    int y = (offset.second + 32) % 64 - 32;
    return std::make_pair(((x + y + 128) / 2) % 64,
                          ((x + y + 128) / 2) % 64);
  }
  default:
    return offset;
  }
}

LifeState Symmetricize(const LifeState &state, StaticSymmetry sym,
                                 std::pair<int, int> offset) {
  switch (sym) {
  case C1:
    return state;
  case C2: {
    LifeState sym = state;
    sym.FlipX();
    sym.FlipY();
    sym.Move(offset.first+1, offset.second+1);
    sym.Join(state);
    return sym;
  }
  case C4: {
    LifeState sym = state;
    sym.Transform(Rotate90);
    sym.Move(offset);
    sym.Join(state);

    LifeState sym2 = sym;
    sym2.FlipX();
    sym2.FlipY();
    sym2.Move(offset.first - offset.second + 1, offset.second + offset.first + 1);
    sym.Join(sym2);
    return sym;
  }
  case D2AcrossX: {
    LifeState sym = state;
    sym.FlipX();
    sym.Move(offset.first, offset.second + 1);
    sym.Join(state);
    return sym;
  }
  case D2AcrossY: {
    LifeState sym = state;
    sym.FlipY();
    sym.Move(offset.first + 1, offset.second);
    sym.Join(state);
    return sym;
  }
  case D2diagodd: {
    LifeState sym = state;
    sym.Transform(ReflectAcrossYeqX);
    sym.Move(offset);
    sym.Join(state);
    return sym;
  }
  case D2negdiagodd: {
    LifeState sym = state;
    sym.Transpose(true);
    sym.Move(offset.first + 1, offset.second + 1);
    sym.Join(state);
    return sym;
  }

  case D4: {
    LifeState acrossx = state;
    auto xoffset = PerpComponent(ReflectAcrossX, offset);
    acrossx.FlipX();
    acrossx.Move(xoffset.first, xoffset.second + 1);
    acrossx.Join(state);

    LifeState acrossy = acrossx;
    auto yoffset = PerpComponent(ReflectAcrossY, offset);
    acrossy.FlipY();
    acrossy.Move(yoffset.first + 1, yoffset.second);
    acrossy.Join(acrossx);

    return acrossy;
  }
  case D4diag: {
    LifeState acrossy = state;
    acrossy.Transform(ReflectAcrossYeqX);
    auto yoffset = PerpComponent(ReflectAcrossYeqX, offset);
    acrossy.Move(yoffset);
    acrossy.Join(state);

    LifeState acrossx = acrossy;
    acrossx.FlipX();
    acrossx.FlipY();
    acrossx.Move(offset.first + 1, offset.second+1);
    acrossx.Join(acrossy);

    return acrossx;
  }

  default:
    __builtin_unreachable();
  }
}

struct SearchState {
  LifeState state;
  LifeState required;
  LifeState history;
  LifeState history1;
  LifeState history2;
  LifeState historyMore;
  LifeState freeState;
  LifeState freeHistory;
  LifeState freeHistory1;
  LifeState freeHistory2;
  LifeState freeHistoryMore;
  StaticSymmetry freeSymmetry;
  std::array<LifeState, 6> triedOffsets;
  Configuration config;

  unsigned endTime;

  unsigned freeCount; // How many catalysts we had at the last free choice
  std::pair<int, int> violationCell;
  unsigned violationGen;

  std::array<int, MAX_CATALYSTS> missingTime;
  std::array<unsigned, MAX_CATALYSTS> recoveredTime;

  std::array<int, MAX_CATALYSTS> freeMissingTime;
  std::array<unsigned, MAX_CATALYSTS> freeRecoveredTime;
};

class CatalystSearcher {
public:
  clock_t begin;
  SearchParams params;
  LifeState pat;
  LifeState alsoRequired;
  std::vector<CatalystData> catalysts;
  std::vector<FilterData> filters;
  std::vector<LifeState> catalystCollisionMasks;

  unsigned nonfixedCatalystCount;
  unsigned fixedCatalystCount;

  unsigned found;
  unsigned fullfound;

  CategoryContainer *categoryContainer;
  CategoryContainer *fullCategoryContainer;

  bool hasFilter;
  bool hasMustInclude;
  bool reportAll;

  unsigned filterMaxGen;

  void LoadMasks() {
    if (params.numCatalysts == 1 || !params.useCollisionMasks)
      return;

    catalystCollisionMasks = std::vector<LifeState>(nonfixedCatalystCount * catalysts.size());

    for (unsigned s = 0; s < catalysts.size(); s++) {
      for (unsigned t = 0; t < nonfixedCatalystCount; t++) {
        if (params.numCatalysts == 2 && hasMustInclude &&
            !catalysts[s].mustInclude && !catalysts[t].mustInclude)
          continue;

        if (params.numTransparent == 1 && catalysts[s].transparent &&
            catalysts[t].transparent)
          continue;

        catalystCollisionMasks[s * nonfixedCatalystCount + t] = CollisionMask(catalysts[s], catalysts[t]);
        catalystCollisionMasks[s * nonfixedCatalystCount + t].RecalculateMinMax();
      }
    }
  }

  void Init(const std::vector<std::string> inputFiles) {
    begin = clock();

    std::vector<CatalystInput> inputcats;
    std::vector<FilterInput> inputfilters;
    for(auto &inputFile : inputFiles)
      ReadParams(inputFile, inputcats, inputfilters, params);

    if (params.pat.length() == 0) {
      std::cout << "Did not read any pattern!" << std::endl;
      exit(1);
    }
    if (inputcats.empty()) {
      std::cout << "Did not read any catalysts!" << std::endl;
      exit(1);
    }

    for (auto &input : inputcats) {
      if(!input.fixed) {
        std::vector<CatalystData> newcats = CatalystData::FromInput(input);
        catalysts.insert(catalysts.end(), newcats.begin(), newcats.end());
      }
    }
    nonfixedCatalystCount = catalysts.size();

    for (auto &input : inputcats) {
      if(input.fixed) {
        std::vector<CatalystData> newcats = CatalystData::FromInput(input);
        catalysts.insert(catalysts.end(), newcats.begin(), newcats.end());
      }
    }
    fixedCatalystCount = catalysts.size() - nonfixedCatalystCount;

    for (auto &input : inputfilters) {
      filters.push_back(FilterData::FromInput(input));
    }

    hasMustInclude = false;
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

    LoadMasks();

    alsoRequired = LifeState::Parse(params.alsoRequired.c_str(), params.alsoRequiredXY.first, params.alsoRequiredXY.second);

    hasFilter = !filters.empty();
    filterMaxGen = FilterMaxGen(filters);

    found = 0;
    fullfound = 0;
    reportAll = params.fullReportFile.length() != 0;
  }

  unsigned FilterMaxGen(std::vector<FilterData> &filters) {
    unsigned maxGen = params.maxGen;

    for (auto &filter : filters) {
      if (filter.gen >= 0 && (unsigned)filter.gen > maxGen)
        maxGen = filter.gen;

      if (filter.range.second >= 0 && (unsigned)filter.range.second > maxGen)
        maxGen = filter.range.second;
    }

    return maxGen;
  }

  void Report(bool saveFile = true) const {
    unsigned sec = (clock() - begin) / CLOCKS_PER_SEC + 1;

    std::cout << "results: " << categoryContainer->categories.size() << "/"
              << found;
    if (params.fullReportFile.length() != 0) {
      std::cout << ", unfiltered: " << fullCategoryContainer->categories.size() << "/"
                << fullfound;
    }
    std::cout << ", now: ";
    PrintTime(sec);
    std::cout << std::endl;

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
    LifeState workspace = Symmetricize(pat, conf.symmetry, conf.symmetryOffset);
    workspace.Join(conf.startingCatalysts);

    std::vector<LifeTarget> movedTargets;

    for (unsigned j = 0; j < params.numCatalysts; j++) {
      for (auto &forbidden : catalysts[conf.curs[j]].forbidden) {
        LifeTarget moved = forbidden;
        moved.wanted.Move(conf.curx[j], conf.cury[j]);
        moved.unwanted.Move(conf.curx[j], conf.cury[j]);
        movedTargets.push_back(moved);
      }
    }

    for (unsigned i = 0; i <= curIter + 1; i++) {
      for (auto &target : movedTargets) {
        if (workspace.Contains(target))
          return false;
      }
      workspace.Step();
    }

    return false;
  }

  bool ValidateFilters(Configuration &conf, unsigned successtime, unsigned failuretime) {
    LifeState workspace = Symmetricize(pat, conf.symmetry, conf.symmetryOffset);
    workspace.Join(conf.startingCatalysts);

    unsigned stopTime;
    if (params.stopAfterCatsDestroyed != -1)
      stopTime = std::min(filterMaxGen, failuretime + params.stopAfterCatsDestroyed);
    else
      stopTime = filterMaxGen;

    LifeState filtersym = Symmetricize(filters[0].target.wanted, conf.symmetry, conf.symmetryOffset);
    filtersym.Step(5);
    LifeTarget filtertarget(filtersym);
    {
      filtersym.Step(5);
      LifeState next = filtersym;
      next.Step(2);
      if(filtersym == next)
        return false;
    }

    for (unsigned g = 0; g <= std::min(filterMaxGen, stopTime); g++) {
      for (auto &filter : filters) {
        bool inSingle = workspace.gen == filter.gen;
        bool inRange = filter.gen == -1 &&
                       filter.range.first <= workspace.gen &&
                       filter.range.second >= workspace.gen;

        bool shouldCheck = inSingle || (inRange && workspace.gen + params.stableInterval >= successtime);

        bool succeeded = false;

        // See whether there is a match at all
        if (shouldCheck && (filter.type == ANDFILTER ||
                            filter.type == ORFILTER)) {
          if(workspace.Contains(filtertarget))
            return true;
        }
      }

      workspace.Step();
    }

    return false;
  }

  bool CheckOscillating(Configuration &conf, unsigned successtime, unsigned failuretime) {
    LifeState workspace = Symmetricize(pat, conf.symmetry, conf.symmetryOffset);
    workspace.Join(conf.startingCatalysts);

    unsigned stopTime;
    if (params.stopAfterCatsDestroyed != -1)
      stopTime = failuretime + params.stopAfterCatsDestroyed;
    else
      stopTime = filterMaxGen;

    std::stack<std::pair<uint64_t, int>> minhashes;

    for(unsigned i = 0; i < std::max((unsigned)50, std::min((unsigned)filters[0].range.second, stopTime)); i++) {
      uint64_t newhash = workspace.GetHash();

      while(true) {
        if(minhashes.empty())
          break;
        if(minhashes.top().first < newhash)
          break;

        if(minhashes.top().first == newhash) {
          unsigned p = i - minhashes.top().second;
          // Avoid some common periods (pentadecathlon and TL hassler)
          return p > 10 && p != 8 && p != 15 && p != 14 && p != 30 && p != 36 && p != 46;
        }

        if(minhashes.top().first > newhash)
          minhashes.pop();
      }

      minhashes.push({newhash, i});

      workspace.Step();
    }

    return false;
  }

  void ReportSolution(Configuration &conf, unsigned successtime, unsigned failuretime) {
    std::pair<int, int> shift = HalveOffset(conf.symmetry, conf.symmetryOffset);
    shift.first = -shift.first;
    shift.second = -shift.second;

    // if reportAll - ignore filters and update fullReport
    if (reportAll) {
      if (HasForbidden(conf, successtime + 3))
        return;

      LifeState workspace = Symmetricize(pat, conf.symmetry, conf.symmetryOffset);
      workspace.Join(conf.startingCatalysts);
      workspace.Move(shift);

      LifeState init = workspace;

      workspace.Step(successtime - params.stableInterval + 2);
      LifeState afterCatalyst = workspace;

      fullfound++;

      fullCategoryContainer->Add(init, afterCatalyst, conf.startingCatalysts, conf,
                                 successtime - params.stableInterval + 2, 0);
    }

    // If has filter validate them;
    if (hasFilter) {
      if (!ValidateFilters(conf, successtime, failuretime)
          && !CheckOscillating(conf, successtime, failuretime))
        return;
    }

    if (HasForbidden(conf, successtime + 3))
      return;

    // If all filters validated update results
    LifeState workspace = Symmetricize(pat, conf.symmetry, conf.symmetryOffset);
    workspace.Join(conf.startingCatalysts);
    workspace.Move(shift);

    LifeState init = workspace;

    workspace.Step(successtime - params.stableInterval + 2);

    LifeState afterCatalyst = workspace;

    categoryContainer->Add(init, afterCatalyst, conf.startingCatalysts, conf,
                           successtime - params.stableInterval + 2, 0);
    found++;
  }

  int OffsetIndexForSym(StaticSymmetry oldsym, StaticSymmetry newsym) {
    if (oldsym != C1) {
      newsym = D2Continuation(oldsym);
    }

    switch (newsym) {
    case C2:
      return 0;
    case C4:
      return 1;
    case D2AcrossX:
      return 2;
    case D2AcrossY:
      return 3;
    case D2diagodd:
      return 4;
    case D2negdiagodd:
      return 5;
    default:
      __builtin_unreachable();
    }
  }

  bool SymIsD2(StaticSymmetry sym) {
    switch (sym){
    case D2AcrossX:
    case D2AcrossY:
    case D2diagodd:
    case D2negdiagodd:
      return true;
    default:
      return false;
    }
  }

  bool SymIsTerminal(StaticSymmetry sym) {
    return sym != C1 && !SymIsD2(sym);
  }

  StaticSymmetry D2Continuation(StaticSymmetry sym) {
    switch (sym){
    case D2AcrossX:
      return D2AcrossY;
    case D2AcrossY:
      return D2AcrossX;
    case D2diagodd:
      return D2negdiagodd;
    case D2negdiagodd:
      return D2diagodd;
    default:
      __builtin_unreachable();
    }
  }

  LifeState IntersectingOffsets(const LifeState &pat1, const LifeState &pat2, StaticSymmetry oldsym, StaticSymmetry newsym) {
    if (oldsym != C1) {
      newsym = D2Continuation(oldsym);
    }

    LifeState transformed;
    switch (newsym) {
    case C2:
      return pat2.Convolve(pat1);
    case C4: {
      transformed = pat1;
      transformed.Transform(Rotate270);

      // // Very inefficient
      // LifeState doubledcollisions = pat2.Convolve(pat1);
      // doubledcollisions.Transform(ReflectAcrossYeqNegXP1);
      // doubledcollisions = doubledcollisions.Skew().HalveY();
      // doubledcollisions.Transform(ReflectAcrossYeqNegXP1);
      // doubledcollisions = doubledcollisions.InvSkew();

      return pat2.Convolve(transformed); // | doubledcollisions;
    }
    case D2AcrossX:
      transformed = pat1;
      transformed.Transform(ReflectAcrossY);
      return pat2.Convolve(transformed);
    case D2AcrossY:
      transformed = pat1;
      transformed.Transform(ReflectAcrossX);
      return pat2.Convolve(transformed);
    case D2diagodd:
      transformed = pat1;
      transformed.Transform(ReflectAcrossYeqNegXP1);
      return pat2.Convolve(transformed);
    case D2negdiagodd:
      transformed = pat1;
      transformed.Transform(ReflectAcrossYeqX);
      return pat2.Convolve(transformed);
    default:
      __builtin_unreachable();
    }
  }

  LifeState IntersectingOffsets(const LifeState &active, StaticSymmetry oldsym,
                                StaticSymmetry newsym) {
    return IntersectingOffsets(active, active, oldsym, newsym);
  }

  LifeState InteractingOffsets(const LifeState &a, const LifeState &a1,
                               const LifeState &a2, const LifeState &aMore,
                               const LifeState &b, const LifeState &b1,
                               const LifeState &b2, const LifeState &bMore,
                               StaticSymmetry oldsym, StaticSymmetry newsym) {
    if (oldsym != C1) {
      newsym = D2Continuation(oldsym);
    }
    switch (newsym) {
    case C2:
    case C4:
      return IntersectingOffsets(a1, b2, oldsym, newsym) |
        IntersectingOffsets(a2, b1, oldsym, newsym) |
        IntersectingOffsets(a | a1 | a2 | aMore, b | bMore, oldsym, newsym) |
        IntersectingOffsets(a | aMore, b | b1 | b2 | bMore, oldsym, newsym);
    case D2AcrossX:
    case D2AcrossY:
    case D2diagodd:
    case D2negdiagodd:
      return IntersectingOffsets(a | a1 | a2 | aMore, b | b1 | b2 | bMore, oldsym, newsym);
    default:
      __builtin_unreachable();
    }
  }

  // Special, faster case when the two patterns are the same
  LifeState InteractingOffsets(const LifeState &a, const LifeState &a1,
                               const LifeState &a2, const LifeState &aMore,
                               StaticSymmetry oldsym, StaticSymmetry newsym) {
    if (oldsym != C1) {
      newsym = D2Continuation(oldsym);
    }
    switch (newsym) {
    case C2:
    case C4:
      return IntersectingOffsets(a1, a2, oldsym, newsym) |
             IntersectingOffsets(a | aMore, a | a1 | a2 | aMore, oldsym, newsym);
    // case C4:
    //   return IntersectingOffsets(a1, a2, oldsym, newsym) |
    //          IntersectingOffsets(a2, a1, oldsym, newsym) |
    //          IntersectingOffsets(a | aMore, a | a1 | a2 | aMore, oldsym, newsym) |
    //          IntersectingOffsets(a | a1 | a2 | aMore, a | aMore, oldsym, newsym);
    case D2AcrossX:
    case D2AcrossY:
    case D2diagodd:
    case D2negdiagodd:
      return IntersectingOffsets(a | a1 | a2 | aMore, a | a1 | a2 | aMore, oldsym, newsym);
    default:
      __builtin_unreachable();
    }
  }

  LifeState AllowedOffsets(StaticSymmetry sym) {
    switch (sym) {
    case C2:
      return ~LifeState();
    case C4:
      return ~LifeState();
    case D2AcrossX:
      return LifeState::SolidRect(0, -32, 1, 64);
    case D2AcrossY:
      return LifeState::SolidRect(-32, 0, 64, 1);
    case D2diagodd: {
      LifeState diag;
      for (int i = 0; i < 64; ++i) {
        diag.SetCell(i, 64 - i, 1);
      }
      return diag;
    }
    case D2negdiagodd: {
      LifeState negdiag;
      for (int i = 0; i < 64; ++i) {
        negdiag.SetCell(i, i, 1);
      }
      return negdiag;
    }
    default:
      __builtin_unreachable();
    }
  }

  std::array<LifeState, 6> StartingOffsets(LifeState &starting) {
    LifeState bounds = LifeState::SolidRect(params.offsetArea[0], params.offsetArea[1], params.offsetArea[2], params.offsetArea[3]);
    std::array<LifeState, 6> result;
    for (auto sym : {C2, C4, D2AcrossX, D2AcrossY, D2diagodd, D2negdiagodd}) {
      result[OffsetIndexForSym(C1, sym)] = IntersectingOffsets(starting, C1, sym) | ~AllowedOffsets(sym) | ~bounds;
    }
    return result;
  }

  void Search() {
    SearchState search;
    search.state = LifeState();
    search.state.JoinWSymChain(pat, params.symmetryChain);
    search.history = search.state;
    search.history1 = search.state;
    search.history2 = search.state;
    search.historyMore = search.state;
    search.freeState = search.state;
    search.freeHistory = search.state;
    search.freeHistory1 = search.state;
    search.freeHistory2 = search.state;
    search.freeHistoryMore = search.state;
    search.freeSymmetry = C1;
    search.freeCount = 0;
    search.freeMissingTime = std::array<int, MAX_CATALYSTS>();
    search.freeMissingTime.fill(-1);
    search.freeRecoveredTime = std::array<unsigned, MAX_CATALYSTS>();
    search.violationCell = {-1, -1};
    search.violationGen = 0;
    search.required = alsoRequired;
    search.missingTime = std::array<int, MAX_CATALYSTS>();
    search.missingTime.fill(-1);
    search.recoveredTime = std::array<unsigned, MAX_CATALYSTS>();
    search.endTime = filterMaxGen;

    LifeState patzoi = pat.ZOI();
    search.triedOffsets = StartingOffsets(patzoi);

    search.config.count = 0;
    search.config.transparentCount = 0;
    search.config.limitedCount = 0;
    search.config.mustIncludeCount = 0;
    search.config.symmetry = C1;

    search.config.startingCatalysts = alsoRequired;

    LifeState bounds =
        LifeState::SolidRect(params.searchArea[0], params.searchArea[1],
                             params.searchArea[2], params.searchArea[3]);

    bounds &= FundamentalDomainFast(params.symmetry);

    std::vector<LifeState> masks = std::vector<LifeState>(nonfixedCatalystCount);
    for (unsigned s = 0; s < nonfixedCatalystCount; s++) {
      LifeState zoi = catalysts[s].state.ZOI();
      zoi.Transform(Rotate180OddBoth);
      masks[s] = search.state.Convolve(zoi) | ~bounds;
    }

    std::array<LifeTarget, MAX_CATALYSTS> shiftedTargets;

    RecursiveSearch(search, masks, shiftedTargets);
  }

  void
  TryApplyingSymmetry(SearchState &search, std::vector<LifeState> &masks,
                      std::array<LifeTarget, MAX_CATALYSTS> &shiftedTargets,
                      const LifeState &activePart,
                      const LifeState &activePart1,
                      const LifeState &activePart2,
                      const LifeState &activePartMore) {

    switch (search.config.symmetry) {
    case C1:
      for (auto sym : {C2, C4, D2AcrossX, D2AcrossY, D2diagodd, D2negdiagodd})
        TryApplyingSpecificSymmetry(search, masks, shiftedTargets, activePart, activePart1, activePart2, activePartMore, sym);
      break;

    case D2AcrossX:
    case D2AcrossY:
      TryApplyingSpecificSymmetry(search, masks, shiftedTargets, activePart, activePart1, activePart2, activePartMore, D4);
      break;

    case D2diagodd:
    case D2negdiagodd:
      TryApplyingSpecificSymmetry(search, masks, shiftedTargets, activePart, activePart1, activePart2, activePartMore, D4diag);
      break;

    default:
      break;
    }
  }

  void TryApplyingSpecificSymmetry(
      SearchState &search, std::vector<LifeState> &masks,
      std::array<LifeTarget, MAX_CATALYSTS> &shiftedTargets,
      const LifeState &state, const LifeState &state1,
      const LifeState &state2, const LifeState &stateMore,
      StaticSymmetry newSym) {

    LifeState newOffsets = InteractingOffsets(state, state1, state2, stateMore,
                                              search.config.symmetry, newSym)
      & ~search.triedOffsets[OffsetIndexForSym(search.config.symmetry, newSym)];
    search.triedOffsets[OffsetIndexForSym(search.config.symmetry, newSym)] |= newOffsets;

    while (!newOffsets.IsEmpty()) {
      // Do the placement
      auto newOffset = newOffsets.FirstOn();

      // HACK: avoid some bad wrapping cases
      if (newSym == D4diag) {
        if (search.config.symmetry == D2diagodd) {
          auto perp = PerpComponent(ReflectAcrossYeqX, newOffset);
          if (!(perp.first == search.config.symmetryOffset.first &&
                perp.second == search.config.symmetryOffset.second)) {
            newOffsets.Erase(newOffset.first, newOffset.second);
            continue;
          }
        }
        if (search.config.symmetry == D2negdiagodd) {
          auto perp = PerpComponent(ReflectAcrossYeqNegXP1, newOffset);
          if (!(perp.first == search.config.symmetryOffset.first &&
                perp.second == search.config.symmetryOffset.second)) {
            newOffsets.Erase(newOffset.first, newOffset.second);
            continue;
          }
        }
      }
      SearchState newSearch = search;
      newSearch.config.symmetry = newSym;
      newSearch.config.symmetryOffset = newOffset;

      newSearch.state = Symmetricize(newSearch.state, newSym, newOffset);

      if (!(search.required & (newSearch.state ^ newSearch.config.startingCatalysts)).IsEmpty()) {
        newOffsets.Erase(newOffset.first, newOffset.second);
        continue;
      }

      {
        LifeState lookahead = newSearch.state;
        lookahead.Step(REQUIRED_LOOKAHEAD);
        if (!(search.required & (lookahead ^ newSearch.config.startingCatalysts)).IsEmpty()) {
          newOffsets.Erase(newOffset.first, newOffset.second);
          continue;
        }
      }

      newSearch.config.startingCatalysts = Symmetricize(newSearch.config.startingCatalysts, newSym, newOffset);
      newSearch.history = Symmetricize(search.history, newSym, newOffset);
      newSearch.history1 = Symmetricize(search.history1, newSym, newOffset);
      newSearch.history2 = Symmetricize(search.history2, newSym, newOffset);
      newSearch.historyMore = Symmetricize(search.historyMore, newSym, newOffset);

      std::vector<LifeState> newMasks;

      if (newSearch.config.count != params.numCatalysts) {
        newMasks = masks;
        LifeState fundamentalDomain = FundamentalDomainFast(newSym);
        fundamentalDomain.Move(HalveOffset(newSym, newOffset));
        for (unsigned t = 0; t < nonfixedCatalystCount; t++) {
          newMasks[t] |= ~fundamentalDomain;
          // At this point, no more catalysts will be placed in this generation.
          // TODO: handle catalysts.reactionMore
          LifeState hitLocations =
              catalysts[t].reaction1.Convolve(newSearch.history2)
            | catalysts[t].reaction2.Convolve(newSearch.history1);

          newMasks[t] |= hitLocations;
        }
      }

      if (search.config.count == 0) {
        std::cout << "Trying offset " << SymmetryToString(newSym) << " "
                  << newOffset.first << ", " << newOffset.second << std::endl;
      }

      if (SymIsD2(newSym)) {
        int continuationIndex = OffsetIndexForSym(C1, D2Continuation(newSym));
        // We need to update the masks for the continuation.
        newSearch.triedOffsets[continuationIndex] = search.triedOffsets[continuationIndex];
        newSearch.triedOffsets[continuationIndex].Move(newOffset);
        // TODO: this is wrong, should use InteractingOffsets somehow:
        newSearch.triedOffsets[continuationIndex] |= IntersectingOffsets(newSearch.history, C1, D2Continuation(newSym));

        if(newSym == D2diagodd || newSym == D2negdiagodd)
          newSearch.triedOffsets[continuationIndex] |= LifeState::Checkerboard();

        RecursiveSearch(newSearch, newMasks, shiftedTargets);
      } else {
        RecursiveSearch(newSearch, newMasks, shiftedTargets);
      }

      newOffsets.Erase(newOffset.first, newOffset.second);
    }
  }

  std::pair<std::pair<int, int>, unsigned> FindViolation(SearchState &search, std::array<LifeTarget, MAX_CATALYSTS> &shiftedTargets) {
    std::array<int, MAX_CATALYSTS> missingTime = search.missingTime;
    LifeState lookahead = search.state;
    for(unsigned g = 1; g <= LIGHTSPEED_LOOKAHEAD; g++) {
      lookahead.Step();

      LifeState requiredViolations = search.required & (lookahead ^ search.config.startingCatalysts);
      std::pair<int, int> cell = requiredViolations.FirstOn();
      if (cell != std::make_pair(-1, -1))
        return {cell, search.state.gen + g};

      for (unsigned i = 0; i < search.config.count; i++) {
        if (lookahead.Contains(shiftedTargets[i]) || catalysts[search.config.curs[i]].sacrificial) {
          missingTime[i] = 0;
        } else {
          missingTime[i] += 1;
        }

        if ((unsigned)missingTime[i] > catalysts[search.config.curs[i]].maxDisappear) {
          std::pair<int, int> cell = (shiftedTargets[i].wanted & ~lookahead).FirstOn();
          if (cell != std::make_pair(-1, -1))
            return {cell, search.state.gen + g};
        }
      }
    }
    return {{-1, -1}, 0};
  }

  void TryAddingFixedCatalyst(SearchState &search, std::vector<LifeState> &masks,
                              std::array<LifeTarget, MAX_CATALYSTS> &shiftedTargets,
                              const LifeState &activeZOI,
                              const LifeState &twonext) {
    for (unsigned s = nonfixedCatalystCount; s < nonfixedCatalystCount + fixedCatalystCount; s++) {
      if (catalysts[s].fixedGen != search.state.gen)
        continue;

      bool alreadyUsed = false;
      for(unsigned i = 0; i < search.config.count; i++) {
        if(search.config.curs[i] == s) {
          alreadyUsed = true;
          break;
        }
      }
      if(alreadyUsed)
        continue;

      if (search.config.transparentCount == params.numTransparent &&
          catalysts[s].transparent)
        continue;
      if (search.config.limitedCount == params.numLimited && catalysts[s].limited)
        continue;
      if (search.config.count == params.numCatalysts - 1 &&
          search.config.mustIncludeCount == 0 && !catalysts[s].mustInclude)
        continue;
      if (!(catalysts[s].state & search.state).IsEmpty())
        continue;
      if((activeZOI & catalysts[s].locusReactionMask).IsEmpty())
        continue;

      // Do the placement
      SearchState newSearch = search;
      newSearch.config.count += 1;
      newSearch.config.curx[search.config.count] = 0;
      newSearch.config.cury[search.config.count] = 0;
      newSearch.config.curs[search.config.count] = s;
      if (catalysts[s].transparent)
        newSearch.config.transparentCount++;
      if (catalysts[s].limited)
        newSearch.config.limitedCount++;
      if (catalysts[s].mustInclude)
        newSearch.config.mustIncludeCount++;

      if(!(catalysts[s].required & search.history).IsEmpty())
        continue;

      newSearch.endTime = std::max(newSearch.endTime, search.state.gen + catalysts[s].maxDisappear);


      LifeState symCatalyst = catalysts[s].state;
      symCatalyst = Symmetricize(symCatalyst, newSearch.config.symmetry,
                                 newSearch.config.symmetryOffset);
      newSearch.config.startingCatalysts |= symCatalyst;
      if(search.state.gen % 2 == 1)
        symCatalyst.Step();
      newSearch.state |= symCatalyst;

      LifeState lookahead = newSearch.state;
      lookahead.Step(2);
      // Do a two-step lookahead to see if the catalyst interacts
      {
        LifeState difference = lookahead ^ (twonext & symCatalyst);
        if (difference.IsEmpty()) {
          continue;
        }
      }

      newSearch.required = search.required | catalysts[s].required;

      for(unsigned i = 0; i < catalysts[s].maxDisappear; i++) {
        lookahead.Step();
        if (!(newSearch.required & (lookahead ^ newSearch.config.startingCatalysts)).IsEmpty()) {
          continue;
        }
      }

      shiftedTargets[search.config.count] = catalysts[s].target;

      if (!lookahead.Contains(shiftedTargets[search.config.count])) {
        continue;
      }

      if (search.config.count == 0 && search.config.symmetry == C1) {
        std::cout << "Placing fixed catalyst " << s << std::endl;
      }

      newSearch.history = search.history | symCatalyst;
      LifeState catalyst1(false), catalyst2(false), catalystMore(false);
      SetCounts(newSearch.config.startingCatalysts, catalyst1, catalyst2, catalystMore);
      // p2 HACK
      if(catalysts[s].periodic) {
        symCatalyst.Step();
        newSearch.history |= symCatalyst;
        UpdateCounts(symCatalyst, catalyst1, catalyst2, catalystMore);
      }

      if(newSearch.config.symmetry == C1) {
        for (auto sym : {C2, C4, D2AcrossX, D2AcrossY, D2diagodd, D2negdiagodd}) {
          unsigned index = OffsetIndexForSym(C1, sym);
          newSearch.triedOffsets[index] |= InteractingOffsets(
              newSearch.history, newSearch.history1, newSearch.history2,
              newSearch.historyMore, symCatalyst, catalyst1, catalyst2,
              catalystMore,
              C1, sym);
        }
      }

      if (SymIsD2(newSearch.config.symmetry)) {
        unsigned index =
            OffsetIndexForSym(newSearch.config.symmetry,
                              D2Continuation(newSearch.config.symmetry));
        newSearch.triedOffsets[index] |= InteractingOffsets(
            newSearch.history, newSearch.history1, newSearch.history2,
            newSearch.historyMore, symCatalyst, catalyst1, catalyst2,
            catalystMore, newSearch.config.symmetry,
            D2Continuation(newSearch.config.symmetry));
      }

      newSearch.history1 |= catalyst1;
      newSearch.history2 |= catalyst2;
      newSearch.historyMore |= catalystMore;

      std::vector<LifeState> newMasks;

      if (newSearch.config.count != params.numCatalysts) {
        UpdateCounts(newSearch.config.startingCatalysts, newSearch.history1, newSearch.history2, newSearch.historyMore);

        newMasks = masks;

        if (params.maxW != -1) {
          LifeState rect =
            LifeState::SolidRect(- params.maxW,
                                 - params.maxH,
                                 2 * params.maxW - 1, 2 * params.maxH - 1);
          LifeState bounds;
          bounds.JoinWSymChain(rect, params.symmetryChain);

          for (unsigned t = 0; t < nonfixedCatalystCount; t++) {
            newMasks[t] |= ~bounds;
          }
        }

        if(params.useCollisionMasks) {
          for (unsigned t = 0; t < nonfixedCatalystCount; t++) {
            newMasks[t].Join(catalystCollisionMasks[s * nonfixedCatalystCount + t]);
          }
        }
      }

      RecursiveSearch(newSearch, newMasks, shiftedTargets);
    }
  }

  void TryAddingCatalyst(SearchState &search, std::vector<LifeState> &masks,
                         std::array<LifeTarget, MAX_CATALYSTS> &shiftedTargets,
                         const LifeState &state1,
                         const LifeState &state2,
                         const LifeState &stateMore,
                         const LifeState &activePart,
                         const LifeState &activePart1,
                         const LifeState &activePart2,
                         const LifeState &activePartMore,
                         const LifeState &twonext,
                         bool justPeriodic) {

    unsigned activePartPop1 = activePart1.GetPop();
    unsigned activePartPop2 = activePart2.GetPop();

    for (unsigned s = 0; s < nonfixedCatalystCount; s++) {
      if (search.config.transparentCount == params.numTransparent &&
          catalysts[s].transparent)
        continue;
      if (search.config.limitedCount == params.numLimited && catalysts[s].limited)
        continue;
      if (search.config.count == params.numCatalysts - 1 &&
          search.config.mustIncludeCount == 0 && !catalysts[s].mustInclude)
        continue;
      if (justPeriodic && !catalysts[s].periodic)
        continue;

      LifeState newPlacements;
      if (!catalysts[s].periodic) {
        if (catalysts[s].hasLocusReactionPop1) {
          if (activePartPop2 < catalysts[s].locusReactionPop1)
            newPlacements |= activePart2.Convolve(catalysts[s].locusReactionMask1);
          else
            newPlacements |= catalysts[s].locusReactionMask1.Convolve(activePart2);
        }
        if (catalysts[s].hasLocusReactionPop2) {
          if (activePartPop1 < catalysts[s].locusReactionPop2)
            newPlacements |= activePart1.Convolve(catalysts[s].locusReactionMask2);
          else
            newPlacements |= catalysts[s].locusReactionMask2.Convolve(activePart1);
        }
        if (catalysts[s].canSmother)
          newPlacements |= activePartMore.Convolve(catalysts[s].locusReactionMask);
      }

      // p2 HACK
      if (catalysts[s].periodic) {
        if (search.state.gen % 2 == catalysts[s].phase){
          newPlacements |= catalysts[s].locusReactionMask1.Convolve(state2);
          newPlacements |= catalysts[s].locusReactionMask2.Convolve(state1);
          if (catalysts[s].canSmother)
            newPlacements |= catalysts[s].locusReactionMask.Convolve(stateMore);
        } else {
          newPlacements |= catalysts[s].offReactionMask1.Convolve(state2);
          newPlacements |= catalysts[s].offReactionMask2.Convolve(state1);
          if (catalysts[s].canSmother)
            newPlacements |= catalysts[s].offReactionMask.Convolve(stateMore);
        }
      }

      newPlacements &= ~masks[s];

      if (!newPlacements.IsEmpty()
          // p2 HACK
          && !catalysts[s].periodic) {
        LifeState hitLocations = catalysts[s].locusAvoidMask.Convolve(activePartMore | search.historyMore);
        hitLocations |= catalysts[s].locusAvoidMask1.Convolve(activePart2 | search.history2);
        hitLocations |= catalysts[s].locusAvoidMask2.Convolve(activePart1 | search.history1);

        masks[s] |= hitLocations;
        newPlacements &= ~hitLocations;
      }

      while (!newPlacements.IsEmpty()) {
        // Do the placement
        auto newPlacement = newPlacements.FirstOn();

        SearchState newSearch = search;
        newSearch.config.count += 1;
        newSearch.config.curx[search.config.count] = newPlacement.first;
        newSearch.config.cury[search.config.count] = newPlacement.second;
        newSearch.config.curs[search.config.count] = s;
        if (catalysts[s].transparent)
          newSearch.config.transparentCount++;
        if (catalysts[s].limited)
          newSearch.config.limitedCount++;
        if (catalysts[s].mustInclude)
          newSearch.config.mustIncludeCount++;

        newSearch.endTime = std::max(newSearch.endTime, search.state.gen + catalysts[s].maxDisappear);

        LifeState shiftedCatalyst = catalysts[s].state;
        shiftedCatalyst.Move(newPlacement.first, newPlacement.second);

        if (!params.useCollisionMasks &&
            !(newSearch.state & shiftedCatalyst).IsEmpty()) {
          masks[s].Set(newPlacement.first, newPlacement.second);
          newPlacements.Erase(newPlacement.first, newPlacement.second);
          continue;
        }

        LifeState symCatalyst = shiftedCatalyst;
        symCatalyst = Symmetricize(symCatalyst, newSearch.config.symmetry,
                                   newSearch.config.symmetryOffset);

        newSearch.state |= symCatalyst;
        newSearch.config.startingCatalysts |= symCatalyst;
        if(catalysts[s].periodic && search.state.gen % 2 == 1) {
            symCatalyst.Step();
        }
        newSearch.state |= symCatalyst;

        LifeState lookahead = newSearch.state;
        lookahead.Step(2);

        // p2 HACK
        // Do a two-step lookahead to see if the catalyst interacts
        {
          LifeState difference = lookahead ^ (twonext & symCatalyst);
          if (difference.IsEmpty()) {
            if (DEBUG_OUTPUT && search.config.count == 0) {
              std::cout << "Skipping catalyst " << s << " at "
                        << newPlacement.first << ", " << newPlacement.second
                        << " (no interaction) " << std::endl;
            }
            // Note: we deliberately don't set the mask,
            // because it may turn out that a catalyst here
            // interacts properly in a later generation.
            newPlacements.Erase(newPlacement.first, newPlacement.second);
            continue;
          }
        }

        if(catalysts[s].hasRequired) {
          newSearch.required.Join(catalysts[s].required, newPlacement.first, newPlacement.second);
        }

        {
          LifeState unused = shiftedCatalyst;
          bool catalystFailed = false;
          for (unsigned i = 2; i < REQUIRED_LOOKAHEAD; i++) {
            unused &= lookahead;

            if ((search.state.gen + i) % 2 == 0 && !(newSearch.required & (lookahead ^ newSearch.config.startingCatalysts)).IsEmpty()) {
              catalystFailed = true;
              break;
            }
            lookahead.Step();
          }
          if (catalystFailed) {
              if (DEBUG_OUTPUT && search.config.count == 0) {
                  std::cout << "Skipping catalyst " << s << " at "
                            << newPlacement.first << ", " << newPlacement.second
                            << " (is destroyed) " << std::endl;
              }
              masks[s].Set(newPlacement.first, newPlacement.second);
              newPlacements.Erase(newPlacement.first, newPlacement.second);
              continue;
          }

          if (!catalysts[s].canRock && unused == shiftedCatalyst) {
            if (DEBUG_OUTPUT && search.config.count == 0) {
              std::cout << "Skipping catalyst " << s << " at "
                        << newPlacement.first << ", " << newPlacement.second
                        << " (is rock) " << std::endl;
            }

            masks[s].Set(newPlacement.first, newPlacement.second);
            newPlacements.Erase(newPlacement.first, newPlacement.second);
            continue;
          }
        }

        shiftedTargets[search.config.count] = catalysts[s].target;
        shiftedTargets[search.config.count].wanted.Move(newPlacement.first, newPlacement.second);
        shiftedTargets[search.config.count].unwanted.Move(newPlacement.first, newPlacement.second);

        if (catalysts[s].checkRecovery) {
          bool catalystFailed = false;
          if(!catalysts[s].checkReaction) {
            int gap = (search.state.gen + catalysts[s].maxDisappear) - lookahead.gen;
            for (int i = 0; i < gap; i++) {
              lookahead.Step();
              // if (catalysts[s].hasRequired && !(catRequired & (lookahead ^ newSearch.config.startingCatalysts)).IsEmpty()) {
              //   catalystFailed = true;
              //   break;
              // }
            }
            if (!lookahead.Contains(shiftedTargets[search.config.count])) {
              catalystFailed = true;
            }
          } else {
            LifeState lookahead = newSearch.state;
            LifeState active;
            for (unsigned i = 0; i < catalysts[s].maxDisappear; i++) {
              lookahead.Step();
              active |= lookahead ^ shiftedCatalyst;
            }
            LifeState unused = catalysts[s].state;
            unused.Move(newPlacement.first, newPlacement.second);
            unused &= ~(newSearch.required | active);
            if (!unused.IsEmpty()) {
              if (DEBUG_OUTPUT && search.config.count == 0) {
                std::cout << "Skipping catalyst " << s << " at "
                          << newPlacement.first << ", " << newPlacement.second
                          << " (failed to react correctly) " << std::endl;
              }

              masks[s].Set(newPlacement.first, newPlacement.second);
              newPlacements.Erase(newPlacement.first, newPlacement.second);
              continue;
            }
            if (!lookahead.Contains(shiftedTargets[search.config.count])) {
              catalystFailed = true;
            }
          }

          if (catalystFailed) {
            if (DEBUG_OUTPUT && search.config.count == 0) {
              std::cout << "Skipping catalyst " << s << " at "
                        << newPlacement.first << ", " << newPlacement.second
                        << " (failed to recover completely) " << std::endl;
            }
            masks[s].Set(newPlacement.first, newPlacement.second);
            newPlacements.Erase(newPlacement.first, newPlacement.second);
            continue;
          }
        }

        if (search.config.count == 0 && search.config.symmetry == C1) {
          std::cout << "Placing catalyst " << s << " at " << newPlacement.first
                    << ", " << newPlacement.second << std::endl;
        }

        LifeState catalyst1(false), catalyst2(false), catalystMore(false);
        if (newSearch.config.count != params.numCatalysts
            || newSearch.config.symmetry == C1
            || SymIsD2(newSearch.config.symmetry)
            ) {
          newSearch.history |= symCatalyst;
          SetCounts(newSearch.config.startingCatalysts, catalyst1, catalyst2, catalystMore);

          // p2 HACK
          if(catalysts[s].periodic) {
            symCatalyst.Step();
            newSearch.history |= symCatalyst;
            UpdateCounts(symCatalyst, catalyst1, catalyst2, catalystMore);
          }

          newSearch.history1 |= catalyst1;
          newSearch.history2 |= catalyst2;
          newSearch.historyMore |= catalystMore;
        }

        if(newSearch.config.symmetry == C1) {
          for (auto sym : {C2, C4, D2AcrossX, D2AcrossY, D2diagodd, D2negdiagodd}) {
            unsigned index = OffsetIndexForSym(C1, sym);
            newSearch.triedOffsets[index] |= InteractingOffsets(
                newSearch.history, newSearch.history1, newSearch.history2,
                newSearch.historyMore, symCatalyst, catalyst1, catalyst2,
                catalystMore, C1, sym);
          }
        }
        if (SymIsD2(newSearch.config.symmetry)) {
          unsigned index =
              OffsetIndexForSym(newSearch.config.symmetry,
                                D2Continuation(newSearch.config.symmetry));
          newSearch.triedOffsets[index] |= InteractingOffsets(
              newSearch.history, newSearch.history1, newSearch.history2,
              newSearch.historyMore, symCatalyst, catalyst1, catalyst2,
              catalystMore, newSearch.config.symmetry,
              D2Continuation(newSearch.config.symmetry));
        }

        std::vector<LifeState> newMasks;

        // If we just placed the last catalyst, don't bother
        // updating the masks
        if (newSearch.config.count != params.numCatalysts) {
          newMasks = masks;
          if (params.maxW != -1) {
            LifeState rect =
                LifeState::SolidRect(newPlacement.first - params.maxW,
                                     newPlacement.second - params.maxH,
                                     2 * params.maxW - 1, 2 * params.maxH - 1);
            LifeState bounds;
            bounds.JoinWSymChain(rect, params.symmetryChain);

            for (unsigned t = 0; t < nonfixedCatalystCount; t++) {
              newMasks[t] |= ~bounds;
            }
          }

          if(params.useCollisionMasks) {
            for (unsigned t = 0; t < nonfixedCatalystCount; t++) {
              newMasks[t].Join(catalystCollisionMasks[s * nonfixedCatalystCount + t],
                               newPlacement.first, newPlacement.second);
            }
          }
        }

        RecursiveSearch(newSearch, newMasks, shiftedTargets);

        masks[s].Set(newPlacement.first, newPlacement.second);
        newPlacements.Erase(newPlacement.first, newPlacement.second);
      }
    }
  }

  void RecursiveSearch(SearchState &search, std::vector<LifeState> &masks, std::array<LifeTarget, MAX_CATALYSTS> &shiftedTargets) {
    bool wasFree = false;
    bool isFree = false;
    if(search.config.count != params.numCatalysts || !SymIsTerminal(search.config.symmetry)) {
      wasFree = search.violationCell.first == -1;
      auto [newViolationCell, newViolationGen] = FindViolation(search, shiftedTargets);
      isFree = newViolationCell.first == -1;

      if (isFree && !wasFree) {
        // Restore free point
        // Note: masks etc are not restored
        LifeState newCatalysts;
        for (unsigned i = search.freeCount; i < search.config.count; i++) {
          LifeState shiftedCatalyst = catalysts[search.config.curs[i]].state;
          if (catalysts[search.config.curs[i]].periodic && search.freeState.gen % 2 == 1)
            shiftedCatalyst.Step();
          shiftedCatalyst.Move(search.config.curx[i], search.config.cury[i]);
          newCatalysts |= shiftedCatalyst;
        }

        if (search.config.symmetry == C1) {
          search.state = search.freeState | newCatalysts;
          search.state.gen = search.freeState.gen;
          search.history = search.freeHistory;

          search.history1 = search.freeHistory1;
          search.history2 = search.freeHistory2;
          search.historyMore = search.freeHistoryMore;
          UpdateCounts(search.config.startingCatalysts, search.history1, search.history2, search.historyMore);
        } else {
          newCatalysts = Symmetricize(newCatalysts, search.config.symmetry, search.config.symmetryOffset);
          search.state = Symmetricize(search.freeState, search.config.symmetry, search.config.symmetryOffset) | newCatalysts;
          search.state.gen = search.freeState.gen;
          search.history = Symmetricize(search.freeHistory, search.config.symmetry, search.config.symmetryOffset);

          search.history1 = Symmetricize(search.freeHistory1, search.config.symmetry, search.config.symmetryOffset);
          search.history2 = Symmetricize(search.freeHistory2, search.config.symmetry, search.config.symmetryOffset);
          search.historyMore = Symmetricize(search.freeHistoryMore, search.config.symmetry, search.config.symmetryOffset);
          UpdateCounts(search.config.startingCatalysts, search.history1, search.history2, search.historyMore);
        }

        search.missingTime = search.freeMissingTime;
        search.recoveredTime = search.freeRecoveredTime;
      }
      search.violationCell = newViolationCell;
      search.violationGen = newViolationGen;
    }

    bool success = false;
    bool failure = false;
    unsigned successtime;
    unsigned failuretime;

    for (unsigned g = search.state.gen; g <= search.endTime; g++) {
      // Block the locations that are hit too early
      if (search.state.gen < params.startGen) {
        search.history |= search.state;

        LifeState state1(false), state2(false), stateMore(false);
        SetCounts(search.state, state1, state2, stateMore);

        search.history1 |= state1;
        search.history2 |= state2;
        search.historyMore |= stateMore;

        for (auto sym :
             {C2, C4, D2AcrossX, D2AcrossY, D2diagodd, D2negdiagodd}) {
          unsigned index = OffsetIndexForSym(C1, sym);
          search.triedOffsets[index] |= InteractingOffsets(
              search.state, state1, state2, stateMore,
              C1, sym);
        }

        search.state.Step();
        continue;
      }
      if (search.state.gen == params.startGen && search.config.count == 0) {
        for (unsigned s = 0; s < nonfixedCatalystCount; s++) {
          masks[s] |= search.history1.Convolve(catalysts[s].locusReactionMask2);
          masks[s] |= search.history2.Convolve(catalysts[s].locusReactionMask1);
          masks[s] |= search.historyMore.Convolve(catalysts[s].locusReactionMask);
        }
      }

      if (search.config.count == 0 && g > params.lastGen)
        failure = true;
      if (search.config.count < params.numCatalysts && g > params.maxGen)
        failure = true;
      if ((search.config.count != params.numCatalysts || !SymIsTerminal(search.config.symmetry)) && !isFree && search.state.gen >= search.violationGen)
        failure = true;
      if (!(search.required & (search.state ^ search.config.startingCatalysts)).IsEmpty())
        failure = true;

      if (failure) {
        failuretime = search.state.gen;
        break;
      }


      if (search.config.count == 0 && search.config.symmetry == C1) {
        std::cout << "Collision at gen " << g << std::endl;
      }

      if(!SymIsTerminal(search.config.symmetry) || search.config.count != params.numCatalysts) {
        LifeState criticalArea(false);
        if (isFree) {
          criticalArea = ~LifeState();
          search.freeState = search.state;
          search.freeHistory = search.history;
          search.freeHistory1 = search.history1;
          search.freeHistory2 = search.history2;
          search.freeHistoryMore = search.historyMore;
          search.freeSymmetry = search.config.symmetry;
          search.freeCount = search.config.count;
          search.freeMissingTime = search.missingTime;
          search.freeRecoveredTime = search.recoveredTime;
        } else {
          int distance = search.violationGen - search.state.gen;
          criticalArea = LifeState::NZOIAround(search.violationCell, distance-1);
          criticalArea = Symmetricize(criticalArea, search.config.symmetry, search.config.symmetryOffset);
        }

        LifeState state1(false), state2(false), stateMore(false);
        SetCounts(search.state, state1, state2, stateMore);

        LifeState newState1 = state1 & ~search.history1 & criticalArea;
        LifeState newState2 = state2 & ~search.history2 & criticalArea;
        LifeState newStateMore = stateMore & ~search.historyMore & criticalArea;

        LifeState allNew = newState1 | newState2 | newStateMore;
        bool hasActivePart = !allNew.IsEmpty();

        // if (hasActivePart) {
          if (search.state.gen >= params.startGen && search.config.count != params.numCatalysts) {
            LifeState next = search.state;
            next.Step();
            LifeState twonext = next;
            twonext.Step();

            TryAddingCatalyst(search, masks, shiftedTargets, state1, state2, stateMore, allNew, newState1, newState2, newStateMore, twonext, !hasActivePart);

            TryAddingFixedCatalyst(search, masks, shiftedTargets, allNew, twonext);
          }
        // }

        TryApplyingSymmetry(search, masks, shiftedTargets, search.state,
                            state1, state2, stateMore);

        search.history |= search.state;
        search.history1 |= state1;
        search.history2 |= state2;
        search.historyMore |= stateMore;

        search.state.Step();
      } else {
        search.state.Step();
      }

      for (unsigned i = 0; i < search.config.count; i++) {
        if (search.recoveredTime[i] == params.stableInterval)
          continue;

        if (search.state.Contains(shiftedTargets[i]) || catalysts[search.config.curs[i]].sacrificial) {
          if (search.missingTime[i] != -1 || catalysts[search.config.curs[i]].canSmother) {
            search.missingTime[i] = 0;
            search.recoveredTime[i] += 1;
          }
        } else {
          // We use -1 as a sentinel to mean the catalyst hasn't interacted yet
          if (search.missingTime[i] == -1)
            search.missingTime[i] = 0;
          search.missingTime[i] += 1;
          search.recoveredTime[i] = 0;
        }

        if (search.missingTime[i] != -1 && (unsigned)search.missingTime[i] > catalysts[search.config.curs[i]].maxDisappear) {
          failuretime = search.state.gen;
          failure = true;
          break;
        }
      }

      if (search.config.count == params.numCatalysts && !success) {
        bool allRecovered = true;
        for (unsigned i = 0; i < search.config.count; i++) {
          if (catalysts[search.config.curs[i]].sacrificial)
            continue;
          if (search.recoveredTime[i] < params.stableInterval || search.missingTime[i] > 0) {
            allRecovered = false;
          }
        }
        if (allRecovered) {
          success = true;
          successtime = g;
        }
      }

      if (search.config.count == 0 && search.config.symmetry == C1)
        Report();
    }

    if (!failure)
      failuretime = filterMaxGen;

    if (success)
      ReportSolution(search.config, successtime, failuretime);
  }
};

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cout << "Usage CatForce.exe <in file>" << std::endl;
    exit(0);
  }

  std::cout << "Input: " << argv[1] << std::endl
            << "Initializing please wait..." << std::endl;

  std::vector<std::string> arguments(argv + 1, argv + argc);

  CatalystSearcher searcher;
  searcher.Init(arguments);

  clock_t initialized = clock();
  printf("Total elapsed time: %f seconds\n",
         (double) (initialized - searcher.begin) / CLOCKS_PER_SEC);
  std::cout << std::endl
            << "Initialization finished, searching..." << std::endl
            << std::endl;

  searcher.Search();

  printf("\n\nFINISH\n");
  clock_t end = clock();
  printf("Total elapsed time: %f seconds\n",
         (double)(end - searcher.begin) / CLOCKS_PER_SEC);
}
