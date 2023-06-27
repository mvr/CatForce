// CatForce - Catalyst search utility based on LifeAPI using brute force.
// Written by Michael Simkin 2015
#include "LifeAPI.h"
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>

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
  int searchArea[4]{};
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
  bool transparent;
  bool limited;
  bool mustInclude;
  bool checkRecovery;
  bool checkReaction;
  bool canSmother;
  bool sacrificial;
  bool fixed;
  int fixedGen;

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
    centerX = atoi(elems[3].c_str());
    centerY = atoi(elems[4].c_str());
    symmType = elems[5].at(0);

    transparent = false;
    limited = false;
    mustInclude = false;
    checkRecovery = false;
    checkReaction = false;
    canSmother = false;
    sacrificial = false;
    fixed = false;
    fixedGen = -1;

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
      } else if (elems[argi] == "sacrificial") {
        sacrificial = true;
        argi += 1;
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

LifeState FundamentalDomain(const StaticSymmetry sym) {
  switch (sym) {
  case StaticSymmetry::C1:
    return LifeState::Parse("64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o!");
  case StaticSymmetry::D2AcrossY:
  case StaticSymmetry::D2AcrossYEven:
    return LifeState::Parse("32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o!");
  case StaticSymmetry::D2AcrossX:
  case StaticSymmetry::D2AcrossXEven:
    return LifeState::Parse("64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o!");
  case StaticSymmetry::D2diagodd:
    return LifeState::Parse("64o$63o$62o$61o$60o$59o$58o$57o$56o$55o$54o$53o$52o$51o$50o$49o$48o$47o$46o$45o$44o$43o$42o$41o$40o$39o$38o$37o$36o$35o$34o$33o$32o$31o$30o$29o$28o$27o$26o$25o$24o$23o$22o$21o$20o$19o$18o$17o$16o$15o$14o$13o$12o$11o$10o$9o$8o$7o$6o$5o$4o$3o$2o$o!");
  case StaticSymmetry::D2negdiagodd:
    return LifeState::Parse("o$2o$3o$4o$5o$6o$7o$8o$9o$10o$11o$12o$13o$14o$15o$16o$17o$18o$19o$20o$21o$22o$23o$24o$25o$26o$27o$28o$29o$30o$31o$32o$33o$34o$35o$36o$37o$38o$39o$40o$41o$42o$43o$44o$45o$46o$47o$48o$49o$50o$51o$52o$53o$54o$55o$56o$57o$58o$59o$60o$61o$62o$63o$64o!");
  case StaticSymmetry::C2:
  case StaticSymmetry::C2even:
  case StaticSymmetry::C2horizontaleven:
  case StaticSymmetry::C2verticaleven:
    return LifeState::Parse("64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o!");
  case StaticSymmetry::C4:
  case StaticSymmetry::C4even:
  case StaticSymmetry::D4:
  case StaticSymmetry::D4even:
  case StaticSymmetry::D4horizontaleven:
  case StaticSymmetry::D4verticaleven:
    return LifeState::Parse("32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o$32o!");
  case StaticSymmetry::D4diag:
  case StaticSymmetry::D4diageven:
    return LifeState::Parse("o$2o$3o$4o$5o$6o$7o$8o$9o$10o$11o$12o$13o$14o$15o$16o$17o$18o$19o$20o$21o$22o$23o$24o$25o$26o$27o$28o$29o$30o$31o$32o$32o$31o$30o$29o$28o$27o$26o$25o$24o$23o$22o$21o$20o$19o$18o$17o$16o$15o$14o$13o$12o$11o$10o$9o$8o$7o$6o$5o$4o$3o$2o$o!");
  case StaticSymmetry::D8:
  case StaticSymmetry::D8even:
    return LifeState::Parse("o$2o$3o$4o$5o$6o$7o$8o$9o$10o$11o$12o$13o$14o$15o$16o$17o$18o$19o$20o$21o$22o$23o$24o$25o$26o$27o$28o$29o$30o$31o$32o!");
  }
}

inline std::pair<int, int> CommuteTranslation(const SymmetryTransform sym, std::pair<int, int> vec) {
  int x = vec.first;
  int y = vec.second;
  switch (sym) {
  case Identity:                return std::make_pair(x  , y);
  case ReflectAcrossXEven:      return std::make_pair(x  , -y);
  case ReflectAcrossX:          return std::make_pair(x  , -y);
  case ReflectAcrossYEven:      return std::make_pair(-x , y);
  case ReflectAcrossY:          return std::make_pair(-x , y);
  case Rotate90Even:            return std::make_pair(-y , x);
  case Rotate90:                return std::make_pair(-y , x);
  case Rotate270Even:           return std::make_pair(y  , -x);
  case Rotate270:               return std::make_pair(y  , -x);
  case Rotate180OddBoth:        return std::make_pair(-x , -y);
  case Rotate180EvenHorizontal: return std::make_pair(-x , -y);
  case Rotate180EvenVertical:   return std::make_pair(-x , -y);
  case Rotate180EvenBoth:       return std::make_pair(-x , -y);
  case ReflectAcrossYeqX:       return std::make_pair(y  , x);
  case ReflectAcrossYeqNegX:    return std::make_pair(-y , -x);
  case ReflectAcrossYeqNegXP1:  return std::make_pair(-y , -x);
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
      std::string symmetryString = "";
      // reverse-compatibility reasons.
      if (elems[1] == "horizontal") {
        symmetryString = "D2|odd";
      } else if (elems[1] == "horizontaleven") {
        symmetryString = "D2|even";
      } else if (elems[1] == "diagonal") {
        symmetryString =
            "D2/"; // I think this was the way that it worked before?
      } else if (elems[1] == "rotate180") {
        symmetryString = "C2";
      } else if (elems[1] == "rotate180evenx") {
        symmetryString = "C2horizontaleven";
      } else if (elems[1] == "rotate180evenboth") {
        symmetryString = "C2evenboth";
      } else {
        symmetryString = elems[1];
      }

      params.symmetry = SymmetryFromString(symmetryString);
      params.symmetryChain = SymmetryChainFromEnum(params.symmetry);
    } else if (elems[0] == alsoRequired) {
      params.alsoRequired = elems[1].c_str();
      params.alsoRequiredXY = std::make_pair(atoi(elems[2].c_str()), atoi(elems[3].c_str()));
    } else if (elems[0] == stopAfterCatsDestroyed){
      params.stopAfterCatsDestroyed = atoi(elems[1].c_str());
    } else if (elems[0] == maxJunk){
      params.maxJunk = atoi(elems[1].c_str());
    } else if (elems[0] == matchSurvive){
      params.matchSurvive = atoi(elems[1].c_str());
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

  bool hasLocusReactionPop;
  bool hasLocusReactionPop1;
  bool hasLocusReactionPop2;
  bool hasLocusReactionPopMore;

  unsigned locusReactionPop;
  unsigned locusReactionPop1;
  unsigned locusReactionPop2;
  unsigned locusReactionPopMore;

  LifeState locusReactionMask;
  LifeState locusReactionMask1;  // Positions that react with a cell with 1 neighbour at the origin etc.
  LifeState locusReactionMask2;
  LifeState locusReactionMaskMore;

  LifeState locusAvoidMask;
  LifeState locusAvoidMask1;
  LifeState locusAvoidMask2;
  LifeState locusAvoidMaskMore;

  LifeState required;

  LifeTarget target;

  LifeState state1;
  LifeState state2;
  LifeState stateMore;

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
  bool sacrificial;
  bool fixed;
  int fixedGen;

  static std::vector<CatalystData> FromInput(CatalystInput &input);
};

std::vector<CatalystData> CatalystData::FromInput(CatalystInput &input) {
  std::vector<SymmetryTransform> trans = CharToTransforms(input.symmType);

  const char *rle = input.rle.c_str();

  std::vector<CatalystData> results;

  for (auto &tran : trans) {
    LifeState pat = LifeState::Parse(rle, input.centerX, input.centerY, tran);

    CatalystData result;

    result.state = pat;
    result.target = LifeTarget(pat);

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

    LifeState locus;
    if (input.locusRLE != "") {
      result.hasLocus = true;
      locus = LifeState::Parse(input.locusRLE.c_str(),
                                      input.locusXY.first,
                                      input.locusXY.second, tran);
    } else {
      result.hasLocus = false;
      locus = pat;
    }

    {
      LifeState shell = result.state.ZOI().Shell();
      LifeState locusZOI = locus.ZOI();
      LifeState nonLocusZOI = (result.state & ~locus).ZOI();

      result.locusAvoidMask = nonLocusZOI & shell;
      result.locusAvoidMask.Transform(Rotate180OddBoth);
      result.locusAvoidMask.RecalculateMinMax();

      result.locusReactionMask = locusZOI & shell;
      result.locusReactionMask.Transform(Rotate180OddBoth);
      result.locusReactionMask.RecalculateMinMax();
      result.locusReactionMask &= ~result.locusAvoidMask;

      // result.locusAvoidMask = result.locusAvoidMask.Shell();

      LifeState state1, state2, stateMore;
      UpdateCounts(result.state, state1, state2, stateMore);

      result.state1 = state1;
      result.state2 = state2;
      result.stateMore = stateMore;

      LifeState badBirths = nonLocusZOI | (~result.state & result.required);

      result.locusReactionMask1 = state1 & shell & ~badBirths;
      result.locusReactionMask1.Transform(Rotate180OddBoth);
      result.locusReactionMask1.RecalculateMinMax();

      result.locusReactionMask2 = state2 & shell & ~badBirths;
      result.locusReactionMask2.Transform(Rotate180OddBoth);
      result.locusReactionMask2.RecalculateMinMax();

      result.locusReactionMaskMore = stateMore & shell & ~nonLocusZOI;
      result.locusReactionMaskMore.Transform(Rotate180OddBoth);
      result.locusReactionMaskMore.RecalculateMinMax();

      result.locusAvoidMask1 = state1 & shell & badBirths;
      result.locusAvoidMask1.Transform(Rotate180OddBoth);
      result.locusAvoidMask1.RecalculateMinMax();

      result.locusAvoidMask2 = state2 & shell & badBirths;
      result.locusAvoidMask2.Transform(Rotate180OddBoth);
      result.locusAvoidMask2.RecalculateMinMax();

      result.locusAvoidMaskMore = stateMore & shell & nonLocusZOI;
      result.locusAvoidMaskMore.Transform(Rotate180OddBoth);
      result.locusAvoidMaskMore.RecalculateMinMax();

      result.locusReactionPop = result.locusReactionMask.GetPop();
      result.locusReactionPop1 = result.locusReactionMask1.GetPop();
      result.locusReactionPop2 = result.locusReactionMask2.GetPop();
      result.locusReactionPopMore = result.locusReactionMaskMore.GetPop();

      result.hasLocusReactionPop = result.locusReactionPop > 0;
      result.hasLocusReactionPop1 = result.locusReactionPop1 > 0;
      result.hasLocusReactionPop2 = result.locusReactionPop2 > 0;
      result.hasLocusReactionPopMore = result.locusReactionPopMore > 0;
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
    result.sacrificial = input.sacrificial;
    result.fixed = input.fixed;
    result.fixedGen = input.fixedGen;

    results.push_back(result);
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
  const char *rle = input.rle.c_str();

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
  std::array<int, MAX_CATALYSTS> curx;
  std::array<int, MAX_CATALYSTS> cury;
  std::array<int, MAX_CATALYSTS> curs;
};

// Fix a, what positions of b causes a collision?
LifeState CollisionMask(const CatalystData &a, const CatalystData &b) {
  LifeState bFlipped = b.state;
  bFlipped.Transform(Rotate180OddBoth);
  LifeState bState1Flipped = b.state1;
  bState1Flipped.Transform(Rotate180OddBoth);
  LifeState bState2Flipped = b.state2;
  bState2Flipped.Transform(Rotate180OddBoth);
  LifeState bStateMoreFlipped = b.stateMore;
  bStateMoreFlipped.Transform(Rotate180OddBoth);

  // Assumes both are still (so will break with periodic)
  LifeState possibleReactionMask =
      (a.state.ZOI().Convolve(bFlipped | bStateMoreFlipped))
    | (a.state1.Convolve(bState2Flipped))
    | (a.state2.Convolve(bState1Flipped));
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


struct SearchState {
  LifeState state;
  LifeState required;
  LifeState history1;
  LifeState history2;
  LifeState historyMore;
  LifeState freeState;
  LifeState freeHistory1;
  LifeState freeHistory2;
  LifeState freeHistoryMore;
  Configuration config;

  unsigned freeCount; // How many catalysts we had at the last free choice
  std::pair<int, int> violationCell;
  unsigned violationGen;

  std::array<unsigned, MAX_CATALYSTS> missingTime;
  std::array<unsigned, MAX_CATALYSTS> recoveredTime;
};

class CatalystSearcher {
public:
  clock_t begin{};
  SearchParams params;
  LifeState pat;
  LifeState alsoRequired;
  std::vector<CatalystData> catalysts;
  std::vector<FilterData> filters;
  std::vector<LifeState> catalystCollisionMasks;

  unsigned nonfixedCatalystCount;
  unsigned fixedCatalystCount;

  unsigned found{};
  unsigned fullfound{};

  CategoryContainer *categoryContainer{};
  CategoryContainer *fullCategoryContainer{};

  bool hasFilter{};
  bool hasMustInclude{};
  bool reportAll{};

  unsigned filterMaxGen{};

  void LoadMasks() {
    if (params.numCatalysts == 1)
      return;

    catalystCollisionMasks = std::vector<LifeState>(catalysts.size() * catalysts.size());

    for (unsigned s = 0; s < catalysts.size(); s++) {
      for (unsigned t = 0; t < catalysts.size(); t++) {
        if (params.numCatalysts == 2 && hasMustInclude &&
            !catalysts[s].mustInclude && !catalysts[t].mustInclude)
          continue;

        if (params.numTransparent == 1 && catalysts[s].transparent &&
            catalysts[t].transparent)
          continue;

        catalystCollisionMasks[s * catalysts.size() + t] = CollisionMask(catalysts[s], catalysts[t]);
        catalystCollisionMasks[s * catalysts.size() + t].RecalculateMinMax();
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
      if (filter.gen >= 0 && filter.gen > maxGen)
        maxGen = filter.gen;

      if (filter.range.second >= 0 && filter.range.second > maxGen)
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
    LifeState workspace = conf.startingCatalysts;
    workspace.JoinWSymChain(pat, params.symmetryChain);

    for (unsigned i = 0; i <= curIter + 1; i++) {
      for (unsigned j = 0; j < params.numCatalysts; j++) {
        for (unsigned k = 0; k < catalysts[conf.curs[j]].forbidden.size(); k++) {
          if (workspace.Contains(catalysts[conf.curs[j]].forbidden[k], conf.curx[j], conf.cury[j]))
            return true;
        }
      }
      workspace.Step();
    }

    return false;
  }

  bool ValidateFilters(Configuration &conf, unsigned successtime, unsigned failuretime) {
    LifeState workspace = conf.startingCatalysts;
    workspace.JoinWSymChain(pat, params.symmetryChain);

    unsigned maxMatchingPop;
    if(params.maxJunk != -1)
      maxMatchingPop = workspace.GetPop() + params.maxJunk;
    else
      maxMatchingPop = 10000;

    std::vector<bool> filterPassed(filters.size(), false);

    unsigned stopTime;
    if (params.stopAfterCatsDestroyed != -1)
      stopTime = std::min(filterMaxGen, failuretime + params.stopAfterCatsDestroyed);
    else
      stopTime = filterMaxGen;

    for (unsigned g = 0; g <= stopTime; g++) {
      unsigned k = 0;
      for (auto &filter : filters) {
        if (filterPassed[k])
          continue; // No need to check it again.

        bool inSingle = workspace.gen == filter.gen;
        bool inRange = filter.gen == -1 &&
                       filter.range.first <= workspace.gen &&
                       filter.range.second >= workspace.gen;
        bool shouldCheck = inSingle || (inRange && workspace.gen + params.stableInterval >= successtime);

        bool succeeded = false;
        LifeState junk;

        // See whether there is a match at all
        if (shouldCheck && (filter.type == ANDFILTER ||
                            filter.type == ORFILTER)) {
          succeeded = workspace.Contains(filter.target);
          if(params.maxJunk != -1)
            junk = workspace & ~filter.target.wanted & ~conf.startingCatalysts;
        }

        if (shouldCheck && (filter.type == MATCHFILTER)) {
          if(workspace.GetPop() <= maxMatchingPop) {
            LifeState withoutCatalysts = workspace & ~conf.startingCatalysts;
            for (auto &target : filter.transformedTargets) {
              LifeState matches = withoutCatalysts.Match(target);
              if (!matches.IsEmpty()) {
                LifeState matchedPart = matches.Convolve(target.wanted);

                // Check that the matched part is not interfered with
                if(params.matchSurvive != -1) {
                  LifeState matchedAdvanced = matchedPart;
                  matchedAdvanced.Step(params.matchSurvive);

                  LifeState workspaceAdvanced = workspace;
                  workspaceAdvanced.Step(params.matchSurvive);

                  if (!(matchedAdvanced & ~workspaceAdvanced).IsEmpty())
                    continue;
                }

                succeeded = true;
                if(params.maxJunk != -1)
                  junk = withoutCatalysts & ~matchedPart;
                break;
              }
            }
          }
        }

        if (succeeded && (params.maxJunk == -1 || junk.GetPop() <= params.maxJunk)) {
          filterPassed[k] = true;

          // If this was an OR filter, consider all the other OR filters passed
          // too.
          if (filter.type == ORFILTER || filter.type == MATCHFILTER) {
            unsigned j = 0;
            for(auto &otherfilter : filters) {
              if (otherfilter.type == ORFILTER || otherfilter.type == MATCHFILTER) {
                filterPassed[j] = true;
              }
              j++;
            }
          }
        }

        // Bail early
        if (workspace.gen == filter.gen &&
            filter.type == ANDFILTER &&
            !workspace.Contains(filter.target)
            )
          return false;
        k++;
      }

      workspace.Step();
    }

    for (unsigned k = 0; k < filters.size(); k++)
      if (!filterPassed[k])
        return false;

    return true;
  }

  void ReportSolution(Configuration &conf, unsigned successtime, unsigned failuretime) {
    if (HasForbidden(conf, successtime + 3))
      return;

    // if reportAll - ignore filters and update fullReport
    if (reportAll) {
      LifeState workspace = conf.startingCatalysts;
      workspace.JoinWSymChain(pat, params.symmetryChain);

      LifeState init = workspace;

      workspace.Step(successtime - params.stableInterval + 2);
      LifeState afterCatalyst = workspace;

      fullfound++;

      fullCategoryContainer->Add(init, afterCatalyst, conf.startingCatalysts, conf,
                                 successtime - params.stableInterval + 2, 0);
    }

    // If has filter validate them;
    if (hasFilter) {
      if (!ValidateFilters(conf, successtime, failuretime))
        return;
    }

    // If all filters validated update results
    LifeState workspace = conf.startingCatalysts;
    workspace.JoinWSymChain(pat, params.symmetryChain);

    LifeState init = workspace;

    workspace.Step(successtime - params.stableInterval + 2);

    LifeState afterCatalyst = workspace;

    categoryContainer->Add(init, afterCatalyst, conf.startingCatalysts, conf,
                           successtime - params.stableInterval + 2, 0);
    found++;
  }

  void Search() {
    SearchState search;
    search.state = LifeState();
    search.state.JoinWSymChain(pat, params.symmetryChain);
    search.history1 = search.state;
    search.history2 = search.state;
    search.historyMore = search.state;
    search.freeState = search.state;
    search.freeHistory1 = search.state;
    search.freeHistory2 = search.state;
    search.freeHistoryMore = search.state;
    search.freeCount = 0;
    search.violationCell = {-1, -1};
    search.violationGen = 0;
    search.required = alsoRequired;
    search.missingTime = std::array<unsigned, MAX_CATALYSTS>();
    search.recoveredTime = std::array<unsigned, MAX_CATALYSTS>();

    search.config.count = 0;
    search.config.transparentCount = 0;
    search.config.limitedCount = 0;
    search.config.mustIncludeCount = 0;

    search.config.startingCatalysts = alsoRequired;

    LifeState bounds =
        LifeState::SolidRect(params.searchArea[0], params.searchArea[1],
                             params.searchArea[2], params.searchArea[3]);

    bounds &= FundamentalDomain(params.symmetry);

    std::vector<LifeState> masks = std::vector<LifeState>(nonfixedCatalystCount);
    for (unsigned s = 0; s < nonfixedCatalystCount; s++) {
      LifeState zoi = catalysts[s].state.ZOI();
      zoi.Transform(Rotate180OddBoth);
      masks[s] = search.state.Convolve(zoi) | ~bounds;
    }

    std::array<LifeTarget, MAX_CATALYSTS> shiftedTargets;

    RecursiveSearch(search, masks, shiftedTargets);
  }

  std::pair<std::pair<int, int>, unsigned> FindViolation(SearchState &search, std::array<LifeTarget, MAX_CATALYSTS> &shiftedTargets) {
    std::array<unsigned, MAX_CATALYSTS> missingTime = search.missingTime;
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

        if (missingTime[i] > catalysts[search.config.curs[i]].maxDisappear) {
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
                              const LifeState &activePart,
                              const LifeState &next) {
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
      if((activePart & catalysts[s].locusReactionMask).IsEmpty())
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

      LifeState symCatalyst = catalysts[s].state;
      symCatalyst.JoinWSymChain(catalysts[s].state, params.symmetryChain);
      newSearch.config.startingCatalysts |= symCatalyst;
      newSearch.state |= symCatalyst;

      LifeState lookahead = newSearch.state;
      lookahead.Step();
      // Do a one-step lookahead to see if the catalyst interacts
      {
        LifeState difference = lookahead ^ next ^ symCatalyst;
        if (difference.IsEmpty()) {
          continue;
        }
      }

      newSearch.required = search.required | catalysts[s].required;

      for(unsigned i = 0; i < REQUIRED_LOOKAHEAD - 1; i++) {
        lookahead.Step();
        if (!(newSearch.required & (lookahead ^ newSearch.config.startingCatalysts)).IsEmpty()) {
          continue;
        }
      }

      shiftedTargets[search.config.count] = catalysts[s].target;

      if (catalysts[s].checkRecovery) {
        bool catalystFailed = false;
        for (int i = 0; i < (int)catalysts[s].maxDisappear - REQUIRED_LOOKAHEAD + 1; i++) {
          lookahead.Step();
          if (!(catalysts[s].required & (lookahead ^ newSearch.config.startingCatalysts)).IsEmpty()) {
            catalystFailed = true;
            break;
          }
        }

        if (!catalystFailed && !lookahead.Contains(shiftedTargets[search.config.count])) {
          catalystFailed = true;
        }

        if (catalystFailed) {
          continue;
        }
      }

      if (search.config.count == 0) {
        std::cout << "Placing fixed catalyst " << s << std::endl;
      }

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

        for (unsigned t = 0; t < nonfixedCatalystCount; t++) {
          newMasks[t].Join(catalystCollisionMasks[s * catalysts.size() + t]);
        }
      }

      RecursiveSearch(newSearch, newMasks, shiftedTargets);
    }
  }

  void TryAddingCatalyst(SearchState &search,
                         std::vector<LifeState> &masks,
                         std::array<LifeTarget, MAX_CATALYSTS> &shiftedTargets,
                         const LifeState &activePart,
                         const LifeState &activePart1,
                         const LifeState &activePart2,
                         const LifeState &activePartMore,
                         const LifeState &next) {

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

      LifeState newPlacements;
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

      newPlacements &= ~masks[s];

      if (!newPlacements.IsEmpty() && catalysts[s].hasLocus) {
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

        LifeState shiftedCatalyst = catalysts[s].state;
        shiftedCatalyst.Move(newPlacement.first, newPlacement.second);

        LifeState symCatalyst;
        symCatalyst.JoinWSymChain(shiftedCatalyst, params.symmetryChain);
        newSearch.state |= symCatalyst;
        newSearch.config.startingCatalysts |= symCatalyst;

        LifeState lookahead = newSearch.state;

        if(catalysts[s].hasRequired) {
          newSearch.required.Join(catalysts[s].required, newPlacement.first, newPlacement.second);
        }

        {
          bool catalystFailed = false;
          const unsigned skip = 1;
          lookahead.Step(skip);
          for (unsigned i = skip; i < REQUIRED_LOOKAHEAD; i++) {
            lookahead.Step();

            if (!(newSearch.required & (lookahead ^ newSearch.config.startingCatalysts)).IsEmpty()) {
              catalystFailed = true;
              break;
            }
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
        }

        shiftedTargets[search.config.count].wanted = shiftedCatalyst;
        shiftedTargets[search.config.count].unwanted = catalysts[s].target.unwanted;
        shiftedTargets[search.config.count].unwanted.Move(newPlacement.first, newPlacement.second);

        if (catalysts[s].checkRecovery) {
          bool catalystFailed = false;
          if(!catalysts[s].checkReaction) {
            int gap = (search.state.gen + catalysts[s].maxDisappear) - (lookahead.gen);
            for (int i = 0; i < gap; i++) {
              lookahead.Step();
              // if (catalysts[s].hasRequired && !(catRequired & (lookahead ^ newSearch.config.startingCatalysts)).IsEmpty()) {
              //   catalystFailed = true;
              //   break;
              // }
            }
          } else {
            // Have to restart for the cells that are used in first lookahead
            lookahead = newSearch.state;
            LifeState unused = shiftedCatalyst & ~newSearch.required;
            for (int i = 0; i < catalysts[s].maxDisappear; i++) {
              lookahead.Step();
              unused &= lookahead;
            }
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
          }

          if (!catalystFailed && !lookahead.Contains(shiftedTargets[search.config.count])) {
            catalystFailed = true;
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

        if (search.config.count == 0) {
          std::cout
              << "Placing catalyst " << s << " at " << newPlacement.first
              << ", " << newPlacement.second << std::endl;
        }

        std::vector<LifeState> newMasks;

        // If we just placed the last catalyst, don't bother
        // updating the masks
        if (newSearch.config.count != params.numCatalysts) {
          UpdateCounts(newSearch.config.startingCatalysts, newSearch.history1, newSearch.history2, newSearch.historyMore);

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

          for (unsigned t = 0; t < nonfixedCatalystCount; t++) {
            newMasks[t].Join(catalystCollisionMasks[s * catalysts.size() + t],
                                    newPlacement.first, newPlacement.second);
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
    if(search.config.count != params.numCatalysts) {
      wasFree = search.violationCell.first == -1;
      auto [newViolationCell, newViolationGen] = FindViolation(search, shiftedTargets);
      isFree = newViolationCell.first == -1;

      if (isFree && !wasFree) {
        // Restore free point
        // Note: masks etc are not restored
        LifeState newCatalysts;
        for (unsigned i = search.freeCount; i < search.config.count; i++) {
          LifeState shiftedCatalyst = catalysts[search.config.curs[i]].state;
          shiftedCatalyst.Move(search.config.curx[i], search.config.cury[i]);
          newCatalysts.JoinWSymChain(shiftedCatalyst, params.symmetryChain);
        }
        search.state = search.freeState | newCatalysts;

        search.history1 = search.freeHistory1;
        search.history2 = search.freeHistory2;
        search.historyMore = search.freeHistoryMore;
        UpdateCounts(search.config.startingCatalysts, search.history1, search.history2, search.historyMore);
      }
      search.violationCell = newViolationCell;
      search.violationGen = newViolationGen;
    }

    bool success = false;
    bool failure = false;
    unsigned successtime;
    unsigned failuretime;

    for (unsigned g = search.state.gen; g <= filterMaxGen; g++) {
      // Block the locations that are hit too early
      if (search.state.gen < params.startGen) {
        UpdateCounts(search.state, search.history1, search.history2, search.historyMore);
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
      if (search.config.count != params.numCatalysts && !isFree && search.state.gen >= search.violationGen)
        failure = true;
      if (!(search.required & (search.state ^ search.config.startingCatalysts)).IsEmpty())
        failure = true;

      if (failure) {
        failuretime = search.state.gen;
        break;
      }

      if (search.config.count == 0) {
        std::cout << "Collision at gen " << g << std::endl;
      }

      if(search.config.count != params.numCatalysts) {
        LifeState criticalArea(false);
        if (isFree) {
          criticalArea = ~LifeState();
          search.freeState = search.state;
          search.freeHistory1 = search.history1;
          search.freeHistory2 = search.history2;
          search.freeHistoryMore = search.historyMore;
          search.freeCount = search.config.count;
        } else {
          int distance = search.violationGen - search.state.gen;
          criticalArea = LifeState::NZOIAround(search.violationCell, distance-1);
          criticalArea.JoinWSymChain(criticalArea, params.symmetryChain);
        }
        LifeState state1, state2, stateMore;
        UpdateCounts(search.state, state1, state2, stateMore);

        LifeState newState1 = state1 & ~search.history1 & criticalArea;
        LifeState newState2 = state2 & ~search.history2 & criticalArea;
        LifeState newStateMore = stateMore & ~search.historyMore & criticalArea;

        LifeState allNew = newState1 | newState2 | newStateMore;
        bool hasActivePart = !allNew.IsEmpty();

        // Try adding a catalyst
        if (hasActivePart && search.state.gen >= params.startGen) {
          LifeState next = search.state;
          next.Step();

          TryAddingCatalyst(search, masks, shiftedTargets, allNew, newState1, newState2, newStateMore, next);

          TryAddingFixedCatalyst(search, masks, shiftedTargets, allNew, next);

          search.history1 |= state1;
          search.history2 |= state2;
          search.historyMore |= stateMore;
          search.state = next;
        } else {
          search.history1 |= state1;
          search.history2 |= state2;
          search.historyMore |= stateMore;
          search.state.Step();
        }
      } else {
        // No need to update history, not used for anything once all
        // catalysts are placed
        search.state.Step();
      }

      for (unsigned i = 0; i < search.config.count; i++) {
        if (search.state.Contains(shiftedTargets[i]) || catalysts[search.config.curs[i]].sacrificial) {
          search.missingTime[i] = 0;
          search.recoveredTime[i] += 1;
        } else {
          search.missingTime[i] += 1;
          search.recoveredTime[i] = 0;
        }

        if (search.missingTime[i] > catalysts[search.config.curs[i]].maxDisappear) {
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

      if (search.config.count == 0)
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
