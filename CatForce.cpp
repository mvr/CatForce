// CatForce - Catalyst search utility based on LifeAPI using brute force.
// Written by Michael Simkin 2015
#include "LifeAPI.h"
#include <algorithm>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <algorithm>

#include <cassert>

const int MAX_CATALYSTS = 5;

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

enum FilterType {
  ANDFILTER,
  ORFILTER,
  MATCHFILTER,
};

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
  std::vector<std::string> targetFilter;
  std::vector<int> filterdx;
  std::vector<int> filterdy;
  std::vector<int> filterGen;
  std::vector<std::pair<int, int>> filterGenRange;
  std::vector<FilterType> filterType;

  int maxCatSize;

  std::string alsoRequired;
  std::pair<int, int> alsoRequiredXY;

  int stopAfterCatsDestroyed;
  int maxJunk;

  int offsetArea[4]{};

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
  bool sacrificial;
  bool periodic;
  bool fixed;

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
    sacrificial = false;
    periodic = false;

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
      } else if (elems[argi] == "sacrificial") {
        sacrificial = true;
        argi += 1;
      } else if (elems[argi] == "period") {
        periodic = true;
        argi += 2;
      } else if (elems[argi] == "fixed") {
        fixed = true;
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

  std::string offsetArea = "offset-area";

  std::string line;

  bool hasLastGen = false;

  while (std::getline(infile, line)) {
    std::vector<std::string> elems = splitwhitespace(line);

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

    if (elems[0] == numLimited)
      params.numLimited = atoi(elems[1].c_str());

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

    if (elems[0] == filter || elems[0] == orfilter || elems[0] == andfilter || elems[0] == matchfilter) {
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

      if (elems[0] == andfilter) {
        params.filterdx.push_back(atoi(elems[3].c_str()));
        params.filterdy.push_back(atoi(elems[4].c_str()));
        params.filterType.push_back(ANDFILTER);
      } else if (elems[0] == orfilter) {
        params.filterdx.push_back(atoi(elems[3].c_str()));
        params.filterdy.push_back(atoi(elems[4].c_str()));
        params.filterType.push_back(ORFILTER);
      } else if (elems[0] == matchfilter) {
        params.filterdx.push_back(0);
        params.filterdy.push_back(0);
        params.filterType.push_back(MATCHFILTER);
      } else {
        // Shouldn't happen
        params.filterdx.push_back(atoi(elems[3].c_str()));
        params.filterdy.push_back(atoi(elems[4].c_str()));
        params.filterType.push_back(ANDFILTER);
      }
    }

    if (elems[0] == maxWH) {
      params.maxW = atoi(elems[1].c_str());
      params.maxH = atoi(elems[2].c_str());
    }

    if (elems[0] == maxCatSize) {
      params.maxCatSize = atoi(elems[1].c_str());
    }

    if (elems[0] == symmetry) {
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
    }
    if (elems[0] == alsoRequired) {
      params.alsoRequired = elems[1].c_str();
      params.alsoRequiredXY = std::make_pair(atoi(elems[2].c_str()), atoi(elems[3].c_str()));
    }
    if (elems[0] == stopAfterCatsDestroyed){
      params.stopAfterCatsDestroyed = atoi(elems[1].c_str());
    }
    if (elems[0] == maxJunk){
      params.maxJunk = atoi(elems[1].c_str());
    }

    if (elems[0] == offsetArea) {
      params.offsetArea[0] = atoi(elems[1].c_str());
      params.offsetArea[1] = atoi(elems[2].c_str());
      params.offsetArea[2] = atoi(elems[3].c_str());
      params.offsetArea[3] = atoi(elems[4].c_str());
    }
  }

  if(!hasLastGen)
    params.lastGen = params.maxGen - 1;
}

class CatalystData {
public:
  LifeState state;
  LifeTarget target;
  LifeState reactionMask;
  unsigned maxDisappear;
  std::vector<LifeTarget> forbidden;
  LifeState required;
  bool hasLocus;
  LifeState locus;
  LifeState locusReactionMask;
  LifeState locusAvoidMask;
  bool transparent;
  bool limited;
  bool mustInclude;
  bool checkRecovery;
  bool sacrificial;
  bool periodic;
  bool fixed;

  static std::vector<CatalystData> FromInput(CatalystInput &input);
};

std::vector<CatalystData> CatalystData::FromInput(CatalystInput &input) {
  std::vector<SymmetryTransform> trans = CharToTransforms(input.symmType);

  const char *rle = input.rle.c_str();

  std::vector<CatalystData> results;

  for (int phase = 0; phase < 2; phase++) {
  for (auto &tran : trans) {
    LifeState pat = LifeState::Parse(rle, input.centerX, input.centerY, tran);
    pat.Step(phase);

    CatalystData result;
    result.state = pat;

    LifeState gen1 = pat;
    gen1.Step();
    result.target = LifeTarget();
    result.target.wanted = pat & gen1;
    result.target.unwanted = (pat | gen1).ZOI() & ~(pat | gen1);

    result.reactionMask = (pat | gen1).BigZOI();
    result.reactionMask.Transform(Rotate180OddBoth);
    result.reactionMask.RecalculateMinMax();

    if (input.locusRLE != "") {
      result.hasLocus = true;
      result.locus = LifeState::Parse(input.locusRLE.c_str(),
                                      input.locusXY.first,
                                      input.locusXY.second, tran);
    } else {
      result.hasLocus = false;
      result.locus = (pat | gen1);
    }

    result.locusReactionMask = result.locus.BigZOI();
    result.locusReactionMask.Transform(Rotate180OddBoth);
    result.locusReactionMask.RecalculateMinMax();

    result.locusAvoidMask = result.reactionMask & ~result.locusReactionMask;
    result.locusAvoidMask.RecalculateMinMax();

    result.maxDisappear = input.maxDisappear;

    for (unsigned k = 0; k < input.forbiddenRLE.size(); k++) {
      result.forbidden.push_back(LifeTarget::Parse(input.forbiddenRLE[k].c_str(),
                                                   input.forbiddenXY[k].first,
                                                   input.forbiddenXY[k].second, tran));
    }

    if (input.requiredRLE != "") {
      result.required |= LifeState::Parse(input.requiredRLE.c_str(),
                                          input.requiredXY.first,
                                          input.requiredXY.second, tran);
    }

    if (input.antirequiredRLE != "") {
      result.required |= LifeState::Parse(input.antirequiredRLE.c_str(),
                                          input.antirequiredXY.first,
                                          input.antirequiredXY.second, tran);
    }

    result.transparent = input.transparent;
    result.limited = input.limited;
    result.mustInclude = input.mustInclude;
    result.checkRecovery = input.checkRecovery;
    result.sacrificial = input.sacrificial;
    result.periodic = input.periodic;
    result.fixed = input.fixed;

    results.push_back(result);
  }
  if(!input.periodic)
    break;
  }
  return results;
}

struct Configuration {
  unsigned count;
  unsigned transparentCount;
  unsigned limitedCount;
  unsigned mustIncludeCount;
  std::array<int, MAX_CATALYSTS> curx;
  std::array<int, MAX_CATALYSTS> cury;
  std::array<int, MAX_CATALYSTS> curs;
  LifeState state;
  LifeState startingCatalysts;

  StaticSymmetry symmetry;
  std::pair<int, int> symmetryOffset;
};

// Fix a, what positions of b causes a collision?
LifeState CollisionMask(const LifeState &a, const LifeState &b) {
  LifeState bReactionMask = b.BigZOI();
  bReactionMask.Transform(Rotate180OddBoth);
  bReactionMask.RecalculateMinMax();
  LifeState possibleReactionMask = bReactionMask.Convolve(a);

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

inline LifeState Symmetricize(const LifeState &state, StaticSymmetry sym,
                              std::pair<int, int> offset) {
  switch (sym) {
  case C1:
    return state;
  case C2: {
    LifeState sym = state;
    sym.Transform(Rotate180OddBoth);
    sym.Move(offset);
    sym.Join(state);
    return sym;
  }
  case C4: {
    LifeState sym = state;
    sym.Transform(Rotate90);
    sym.Move(offset);
    sym.Join(state);
    sym.Transform(Rotate90);
    sym.Move(offset);
    sym.Join(state);
    sym.Transform(Rotate90);
    sym.Move(offset);
    sym.Join(state);
    return sym;
  }
  case D2AcrossX: {
    LifeState sym = state;
    sym.Transform(ReflectAcrossX);
    sym.Move(offset);
    sym.Join(state);
    return sym;
  }
  case D2AcrossY: {
    LifeState sym = state;
    sym.Transform(ReflectAcrossY);
    sym.Move(offset);
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
    sym.Transform(ReflectAcrossYeqNegXP1);
    sym.Move(offset);
    sym.Join(state);
    return sym;
  }
  case D4: {
    LifeState acrossx = state;
    acrossx.Transform(ReflectAcrossX);
    auto xoffset = PerpComponent(ReflectAcrossX, offset);
    acrossx.Move(xoffset);
    acrossx.Join(state);
    LifeState acrossy = acrossx;
    acrossy.Transform(ReflectAcrossY);
    auto yoffset = PerpComponent(ReflectAcrossY, offset);
    acrossy.Move(yoffset);
    acrossy.Join(acrossx);
    return acrossy;
  }
  case D4diag: {
    LifeState acrossx = state;
    acrossx.Transform(ReflectAcrossYeqNegXP1);
    auto xoffset = PerpComponent(ReflectAcrossYeqNegXP1, offset);
    acrossx.Move(xoffset);
    acrossx.Join(state);
    LifeState acrossy = acrossx;
    acrossy.Transform(ReflectAcrossYeqX);
    auto yoffset = PerpComponent(ReflectAcrossYeqX, offset);
    acrossy.Move(yoffset);
    acrossy.Join(acrossx);
    return acrossy;
  }

  default:
    __builtin_unreachable();
  }
}

class CatalystSearcher {
public:
  clock_t begin{};
  SearchParams params;
  LifeState pat;
  LifeState alsoRequired;
  std::vector<CatalystData> catalysts;
  std::vector<LifeTarget> targetFilter;
  std::vector<LifeState> catalystCollisionMasks;

  unsigned found{};
  unsigned fullfound{};

  CategoryContainer *categoryContainer{};
  CategoryContainer *fullCategoryContainer{};

  bool hasFilter{};
  bool hasMustInclude{};
  bool reportAll{};

  unsigned filterMaxGen{};

  uint64_t AllCatalystsHash() const {
    uint64_t result = 0;
    for (auto &cat : catalysts) {
      result ^= cat.state.GetHash();
    }
    return result;
  }

  void LoadMasks() {
    if (params.numCatalysts == 1)
      return;

    catalystCollisionMasks = std::vector<LifeState>(catalysts.size() * catalysts.size());

    std::stringstream ss;
    ss << "maskpack-" << AllCatalystsHash();
    std::string fname = ss.str();

    // See if the pack exists
    std::ifstream infile;
    infile.open(fname.c_str(), std::ios::binary);
    if (infile.good()) {
      for (unsigned s = 0; s < catalysts.size(); s++) {
        for (unsigned t = 0; t < catalysts.size(); t++) {
          infile.read((char*)catalystCollisionMasks[s * catalysts.size() + t].state, N * sizeof(uint64_t));
          catalystCollisionMasks[s * catalysts.size() + t].RecalculateMinMax();
        }
      }
      return;
    }

    // If not, load or generate the masks
    for (unsigned s = 0; s < catalysts.size(); s++) {
      for (unsigned t = 0; t < catalysts.size(); t++) {
        if (params.numCatalysts == 2 && hasMustInclude &&
            !catalysts[s].mustInclude && !catalysts[t].mustInclude)
          continue;

        if (params.numTransparent == 1 && catalysts[s].transparent &&
            catalysts[t].transparent)
          continue;

        if(catalysts[s].sacrificial || catalysts[t].sacrificial)
          continue;

        catalystCollisionMasks[s * catalysts.size() + t] = CollisionMask(catalysts[s].state, catalysts[t].state);
        catalystCollisionMasks[s * catalysts.size() + t].RecalculateMinMax();
      }
    }

    // Save the pack for next time
    std::ofstream outfile;
    outfile.open(fname.c_str(), std::ofstream::binary);
    for (unsigned s = 0; s < catalysts.size(); s++) {
      for (unsigned t = 0; t < catalysts.size(); t++) {
        outfile.write((char *)catalystCollisionMasks[s * catalysts.size() + t].state, N * sizeof(uint64_t));
      }
    }
    outfile.close();
  }

  void Init(const std::vector<std::string> inputFiles) {
    begin = clock();

    std::vector<CatalystInput> inputcats;
    for(auto &inputFile : inputFiles)
      ReadParams(inputFile, inputcats, params);

    if (params.pat.length() == 0) {
      std::cout << "Did not read any pattern!" << std::endl;
      exit(1);
    }
    if (inputcats.empty()) {
      std::cout << "Did not read any catalysts!" << std::endl;
      exit(1);
    }

    for (auto &input : inputcats) {
      std::vector<CatalystData> newcats = CatalystData::FromInput(input);
      catalysts.insert(catalysts.end(), newcats.begin(), newcats.end());
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

    for (unsigned i = 0; i < params.targetFilter.size(); i++)
      targetFilter.push_back(LifeTarget::Parse(params.targetFilter[i].c_str(),
                                               params.filterdx[i], params.filterdy[i]));

    LoadMasks();

    alsoRequired = LifeState::Parse(params.alsoRequired.c_str(), params.alsoRequiredXY.first, params.alsoRequiredXY.second);

    found = 0;
    fullfound = 0;

    hasFilter = !params.targetFilter.empty();
    reportAll = params.fullReportFile.length() != 0;

    filterMaxGen = FilterMaxGen();
  }

  unsigned FilterMaxGen() {
    unsigned maxGen = params.maxGen;

    for (unsigned j = 0; j < targetFilter.size(); j++) {
      if (params.filterGen[j] >= 0 && params.filterGen[j] > maxGen)
        maxGen = params.filterGen[j];

      if (params.filterGenRange[j].second >= 0 && params.filterGenRange[j].second > maxGen)
        maxGen = params.filterGenRange[j].second;
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

  bool FilterForCurrentGenFail(LifeState &workspace) {
    for (unsigned j = 0; j < targetFilter.size(); j++) {
      if (workspace.gen == params.filterGen[j] &&
          workspace.Contains(targetFilter[j]) == false) {
        return true;
      }
    }

    return false;
  }

  bool ValidateFilters(Configuration &conf, unsigned successtime, unsigned failuretime) {
    LifeState workspace = Symmetricize(pat, conf.symmetry, conf.symmetryOffset);
    workspace.Join(conf.startingCatalysts);

    unsigned stopTime;
    if (params.stopAfterCatsDestroyed != -1)
      stopTime = failuretime + params.stopAfterCatsDestroyed;
    else
      stopTime = filterMaxGen;

    LifeState filtersym = Symmetricize(targetFilter[0].wanted, conf.symmetry, conf.symmetryOffset);
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
      for (unsigned k = 0; k < params.filterGen.size(); k++) {
        bool inSingle = workspace.gen == params.filterGen[k];
        bool inRange = params.filterGen[k] == -1 &&
                       params.filterGenRange[k].first <= workspace.gen &&
                       params.filterGenRange[k].second >= workspace.gen;
        bool shouldCheck = inSingle || (inRange && workspace.gen + params.stableInterval >= successtime);

        bool succeeded = false;

        // See whether there is a match at all
        if (shouldCheck && (params.filterType[k] == ANDFILTER ||
                            params.filterType[k] == ORFILTER)) {
          if(workspace.Contains(filtertarget))
            return true;
        }
      }

      workspace.Step();
    }

    return false;
  }

  void ReportSolution(Configuration &conf, unsigned successtime, unsigned failuretime) {
    if (HasForbidden(conf, successtime + 3))
      return;

    std::pair<int, int> shift = HalveOffset(conf.symmetry, conf.symmetryOffset);
    shift.first = -shift.first;
    shift.second = -shift.second;

    // if reportAll - ignore filters and update fullReport
    if (reportAll) {
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
      if (!ValidateFilters(conf, successtime, failuretime))
        return;
    }

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
      exit(1);
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
      exit(1);
      __builtin_unreachable();
    }
  }

  LifeState CollidingOffsets(LifeState &active, StaticSymmetry oldsym,
                       StaticSymmetry newsym) {
    if (oldsym != C1) {
      newsym = D2Continuation(oldsym);
    }

    LifeState transformed;
    switch (newsym) {
    case C2:
      return active.Convolve(active);
    case C4: {
      transformed = active;
      transformed.Transform(Rotate270);

      // Very inefficient
      LifeState doubledcollisions = active.Convolve(active);
      doubledcollisions.Transform(ReflectAcrossYeqNegXP1);
      doubledcollisions = doubledcollisions.Skew().HalveY();
      doubledcollisions.Transform(ReflectAcrossYeqNegXP1);
      doubledcollisions = doubledcollisions.InvSkew();

      return active.Convolve(transformed) | doubledcollisions;
    }
    case D2AcrossX:
      transformed = active;
      transformed.Transform(ReflectAcrossY);
      return active.Convolve(transformed);
    case D2AcrossY:
      transformed = active;
      transformed.Transform(ReflectAcrossX);
      return active.Convolve(transformed);
    case D2diagodd:
      transformed = active;
      transformed.Transform(ReflectAcrossYeqNegXP1);
      return active.Convolve(transformed);
    case D2negdiagodd:
      transformed = active;
      transformed.Transform(ReflectAcrossYeqX);
      return active.Convolve(transformed);
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
      exit(1);
      __builtin_unreachable();
    }
  }

  std::array<LifeState, 6> StartingOffsets(LifeState &starting) {
    LifeState bounds = LifeState::SolidRect(params.offsetArea[0], params.offsetArea[1], params.offsetArea[2], params.offsetArea[3]);
    std::array<LifeState, 6> result;
    for (auto sym : {C2, C4, D2AcrossX, D2AcrossY, D2diagodd, D2negdiagodd}) {
      result[OffsetIndexForSym(C1, sym)] = CollidingOffsets(starting, C1, sym) | ~AllowedOffsets(sym) | ~bounds;
    }
    return result;
  }

  void Search() {
    Configuration config;
    config.count = 0;
    config.transparentCount = 0;
    config.limitedCount = 0;
    config.mustIncludeCount = 0;
    config.state.JoinWSymChain(pat, params.symmetryChain);
    config.symmetry = C1;

    LifeState history = config.state;

    LifeState bounds =
        LifeState::SolidRect(params.searchArea[0], params.searchArea[1],
                             params.searchArea[2], params.searchArea[3]);

    bounds &= FundamentalDomainFast(params.symmetry);

    LifeState pinhole = ~LifeState();
    pinhole.Erase(0, 0);

    std::vector<LifeState> masks(catalysts.size());
    for (unsigned s = 0; s < catalysts.size(); s++) {
      LifeState zoi = catalysts[s].state.ZOI();
      zoi.Transform(Rotate180OddBoth);
      masks[s] = config.state.Convolve(zoi) | ~bounds;
      if(catalysts[s].fixed) {
        masks[s] = pinhole;
      }
    }

    LifeState patzoi = pat.ZOI();
    std::array<LifeState, 6> triedOffsets = StartingOffsets(patzoi);

    std::vector<LifeTarget> shiftedTargets(params.numCatalysts);

    RecursiveSearch(config, history, alsoRequired, masks,
                    shiftedTargets, triedOffsets,
                    std::array<unsigned, MAX_CATALYSTS>(),
                    std::array<unsigned, MAX_CATALYSTS>());
  }

  void TryApplyingSymmetry(
      Configuration &config, LifeState &history, const LifeState &required,
      std::vector<LifeState> &masks,
      std::vector<LifeTarget> &shiftedTargets, // This can be shared

      std::array<LifeState, 6> &triedOffsets,

      std::array<unsigned, MAX_CATALYSTS> &missingTime,
      std::array<unsigned, MAX_CATALYSTS> &recoveredTime,

      LifeState &activePart) {

    switch (config.symmetry) {
    case C1:
      for (auto sym : {C2, C4, D2AcrossX, D2AcrossY, D2diagodd, D2negdiagodd})
        TryApplyingSpecificSymmetry(
            config, history, required, masks, shiftedTargets,
            triedOffsets, missingTime, recoveredTime, activePart, sym);
      break;

    case D2AcrossX:
    case D2AcrossY:
      TryApplyingSpecificSymmetry(
          config, history, required, masks, shiftedTargets,
          triedOffsets, missingTime, recoveredTime, activePart, D4);
      break;

    case D2diagodd:
    case D2negdiagodd:
      TryApplyingSpecificSymmetry(
          config, history, required, masks, shiftedTargets,
          triedOffsets, missingTime, recoveredTime, activePart, D4diag);
      break;

    default:
      break;
    }
  }

  void TryApplyingSpecificSymmetry(
      Configuration &config, LifeState &history, const LifeState &required,
      std::vector<LifeState> &masks,
      std::vector<LifeTarget> &shiftedTargets, // This can be shared

      std::array<LifeState, 6> &triedOffsets,

      std::array<unsigned, MAX_CATALYSTS> &missingTime,
      std::array<unsigned, MAX_CATALYSTS> &recoveredTime,

      LifeState &activePart, StaticSymmetry newSym) {

    LifeState activezoi = activePart.ZOI();
    LifeState newOffsets = CollidingOffsets(activezoi, config.symmetry, newSym) &
      ~triedOffsets[OffsetIndexForSym(config.symmetry, newSym)];
    triedOffsets[OffsetIndexForSym(config.symmetry, newSym)] |= newOffsets;

    while (!newOffsets.IsEmpty()) {
      // Do the placement
      auto newOffset = newOffsets.FirstOn();

      // HACK: avoid some bad wrapping cases
      if (newSym == D4diag) {
        if (config.symmetry == D2diagodd) {
          auto perp = PerpComponent(ReflectAcrossYeqX, newOffset);
          if (!(perp.first == config.symmetryOffset.first &&
                perp.second == config.symmetryOffset.second)) {
            newOffsets.Erase(newOffset.first, newOffset.second);
            continue;
          }
        }
        if (config.symmetry == D2negdiagodd) {
          auto perp = PerpComponent(ReflectAcrossYeqNegXP1, newOffset);
          if (!(perp.first == config.symmetryOffset.first &&
                perp.second == config.symmetryOffset.second)) {
            newOffsets.Erase(newOffset.first, newOffset.second);
            continue;
          }
        }
      }

      Configuration newConfig = config;
      newConfig.symmetry = newSym;
      newConfig.symmetryOffset = newOffset;

      newConfig.state = Symmetricize(newConfig.state, newSym, newOffset);

      {
        LifeState lookahead = config.state;
        lookahead.Step(10);
        if (!(required & (lookahead ^ newConfig.startingCatalysts)).IsEmpty()) {
          newOffsets.Erase(newOffset.first, newOffset.second);
          continue;
        }
      }

      newConfig.startingCatalysts =
          Symmetricize(newConfig.startingCatalysts, newSym, newOffset);

      LifeState newHistory = Symmetricize(history, newSym, newOffset);

      std::vector<LifeState> newMasks = masks;

      if (newConfig.count != params.numCatalysts) {
        LifeState fundamentalDomain = FundamentalDomainFast(newSym);
        fundamentalDomain.Move(HalveOffset(newSym, newOffset));
        for (unsigned t = 0; t < catalysts.size(); t++) {
          if(!catalysts[t].fixed)
            newMasks[t] |= ~fundamentalDomain | newHistory.Convolve(catalysts[t].reactionMask);
        }
      }

      if (config.count == 0) {
        std::cout << "Trying offset " << SymmetryToString(newSym) << " "
                  << newOffset.first << ", " << newOffset.second << std::endl;
      }

      if (SymIsD2(newSym)) {
        int continuationIndex = OffsetIndexForSym(C1, D2Continuation(newSym));
        // We need to update the masks for the continuation.
        std::array<LifeState, 6> newTriedOffsets;
        newTriedOffsets[continuationIndex] = triedOffsets[continuationIndex];
        newTriedOffsets[continuationIndex].Move(newOffset);
        newTriedOffsets[continuationIndex] |= CollidingOffsets(newHistory, C1, D2Continuation(newSym));

        if(newSym == D2diagodd || newSym == D2negdiagodd)
          newTriedOffsets[continuationIndex] |= LifeState::Checkerboard();

        RecursiveSearch(newConfig, newHistory, required, newMasks,
                        shiftedTargets, newTriedOffsets, missingTime, recoveredTime);
      } else {
        RecursiveSearch(newConfig, newHistory, required, newMasks,
                        shiftedTargets, triedOffsets, missingTime, recoveredTime);
      }

      newOffsets.Erase(newOffset.first, newOffset.second);
    }
  }

  void TryAddingCatalyst(
      Configuration &config, LifeState &history, const LifeState &required,
      std::vector<LifeState> &masks,
      std::vector<LifeTarget> &shiftedTargets, // This can be shared

      std::array<LifeState, 6> &triedOffsets,

      std::array<unsigned, MAX_CATALYSTS> &missingTime,
      std::array<unsigned, MAX_CATALYSTS> &recoveredTime, LifeState &activePart,
      LifeState &twonext) {

    for (unsigned s = 0; s < catalysts.size(); s++) {
      if (config.transparentCount == params.numTransparent &&
          catalysts[s].transparent)
        continue;
      if (config.limitedCount == params.numLimited && catalysts[s].limited)
        continue;
      if (config.count == params.numCatalysts - 1 &&
          config.mustIncludeCount == 0 && !catalysts[s].mustInclude)
        continue;

      LifeState newPlacements =
          activePart.Convolve(catalysts[s].locusReactionMask) & ~masks[s];

      while (!newPlacements.IsEmpty()) {
        // Do the placement
        auto newPlacement = newPlacements.FirstOn();

        Configuration newConfig = config;
        newConfig.count += 1;
        newConfig.curx[config.count] = newPlacement.first;
        newConfig.cury[config.count] = newPlacement.second;
        newConfig.curs[config.count] = s;
        if (catalysts[s].transparent)
          newConfig.transparentCount++;
        if (catalysts[s].limited)
          newConfig.limitedCount++;
        if (catalysts[s].mustInclude)
          newConfig.mustIncludeCount++;

        LifeState shiftedCatalyst = catalysts[s].state;
        shiftedCatalyst.Move(newPlacement.first, newPlacement.second);

        LifeState symCatalyst = shiftedCatalyst;
        symCatalyst = Symmetricize(symCatalyst, newConfig.symmetry,
                                   newConfig.symmetryOffset);
        newConfig.startingCatalysts |= symCatalyst;
        if(config.state.gen % 2 == 1)
          symCatalyst.Step();
        newConfig.state |= symCatalyst;

        // Do a two-step lookahead to see if the catalyst interacts
        {
          LifeState newnext = newConfig.state;
          newnext.Step(2);

          LifeState difference = newnext ^ twonext ^ symCatalyst;
          if (difference.IsEmpty()) {
            // if (config.count == 0 && config.symmetry == C1) {
            //   std::cout << "Skipping catalyst " << s << " at "
            //             << newPlacement.first << ", " << newPlacement.second
            //             << " (no interaction) " << std::endl;
            // }

            // Note: we deliberately don't set the mask,
            // because it may turn out that a catalyst here
            // interacts properly in a later generation.
            newPlacements.Erase(newPlacement.first, newPlacement.second);
            continue;
          }
        }

        LifeState newRequired = required;
        newRequired.Join(catalysts[s].required, newPlacement.first,
                         newPlacement.second);

        {
          LifeState lookahead = newConfig.state;
          lookahead.Step(10);
          if (!(newRequired & (lookahead ^ newConfig.startingCatalysts)).IsEmpty()) {
            // if (config.count == 0 && config.symmetry == C1) {
            //   std::cout << "Skipping catalyst " << s << " at "
            //             << newPlacement.first << ", " << newPlacement.second
            //             << " (is destroyed) " << std::endl;
            // }

            masks[s].Set(newPlacement.first, newPlacement.second);
            newPlacements.Erase(newPlacement.first, newPlacement.second);
            continue;
          }
        }

        shiftedTargets[config.count] = catalysts[s].target;
        shiftedTargets[config.count].wanted.Move(newPlacement.first, newPlacement.second);
        shiftedTargets[config.count].unwanted.Move(newPlacement.first, newPlacement.second);

        if (catalysts[s].checkRecovery) {
          LifeState lookahead = newConfig.state;
          lookahead.Step(catalysts[s].maxDisappear);

          if (!lookahead.Contains(shiftedTargets[config.count])) {
            // if (config.count == 0 && config.symmetry == C1) {
            //   std::cout << "Skipping catalyst " << s << " at "
            //             << newPlacement.first << ", " << newPlacement.second
            //             << " (failed to recover completely) " << std::endl;
            // }

            masks[s].Set(newPlacement.first, newPlacement.second);
            newPlacements.Erase(newPlacement.first, newPlacement.second);
            continue;
          }
        }

        if (config.count == 0 && config.symmetry == C1) {
          std::cout << "Placing catalyst " << s << " at " << newPlacement.first
                    << ", " << newPlacement.second << std::endl;
        }

        LifeState newHistory = history | symCatalyst;
        symCatalyst.Step();
        newHistory |= symCatalyst;

        std::array<LifeState, 6> newOffsets = triedOffsets;
        if(newConfig.symmetry == C1) {
          for (auto sym : {C2, C4, D2AcrossX, D2AcrossY, D2diagodd, D2negdiagodd}) {
            newOffsets[OffsetIndexForSym(C1, sym)] |= CollidingOffsets(newHistory, C1, sym);
          }
        }
        if (SymIsD2(newConfig.symmetry)) {
          newOffsets[OffsetIndexForSym(newConfig.symmetry, D2Continuation(newConfig.symmetry))] |= CollidingOffsets(newHistory, newConfig.symmetry, D2Continuation(newConfig.symmetry));
        }

        std::vector<LifeState> newMasks = masks;

        // If we just placed the last catalyst, don't bother
        // updating the masks
        if (newConfig.count != params.numCatalysts) {
          LifeState bounds;
          if (params.maxW != -1) {
            LifeState rect =
                LifeState::SolidRect(newPlacement.first - params.maxW,
                                     newPlacement.second - params.maxH,
                                     2 * params.maxW - 1, 2 * params.maxH - 1);
            bounds.JoinWSymChain(rect, params.symmetryChain);

            for (unsigned t = 0; t < catalysts.size(); t++) {
              newMasks[t] |= ~bounds;
            }
          }

          for (unsigned t = 0; t < catalysts.size(); t++) {
            newMasks[t].Join(catalystCollisionMasks[s * catalysts.size() + t],
                             newPlacement.first, newPlacement.second);
          }
        }

        RecursiveSearch(newConfig, newHistory, newRequired,
                        newMasks, shiftedTargets, newOffsets, missingTime,
                        recoveredTime);

        masks[s].Set(newPlacement.first, newPlacement.second);
        newPlacements.Erase(newPlacement.first, newPlacement.second);
      }
    }

    history |= config.state;
  }

  void
  RecursiveSearch(Configuration config, LifeState history, const LifeState required,
                  std::vector<LifeState> masks,
                  std::vector<LifeTarget> &shiftedTargets, // This can be shared

                  std::array<LifeState, 6> &triedOffsets,

                  std::array<unsigned, MAX_CATALYSTS> missingTime,
                  std::array<unsigned, MAX_CATALYSTS> recoveredTime) {
    bool success = false;
    bool failure = false;
    unsigned successtime;
    unsigned failuretime;

    for (unsigned g = config.state.gen; g < filterMaxGen; g++) {
      // Block the locations that are hit too early
      if (config.state.gen < params.startGen) {
        for (unsigned s = 0; s < catalysts.size(); s++) {
          LifeState hitLocations = config.state.Convolve(catalysts[s].reactionMask);
          masks[s] |= hitLocations;
        }
        config.state.Step();
        continue;
      }

      if (config.count == 0 && g > params.lastGen)
        failure = true;
      if (config.count < params.numCatalysts && g > params.maxGen)
        failure = true;
      if (!(required & (config.state ^ config.startingCatalysts)).IsEmpty())
        failure = true;

      if (failure) {
        failuretime = config.state.gen;
        break;
      }

      if (config.count == 0 && config.symmetry == C1) {
        std::cout << "Collision at gen " << g << std::endl;
      }

      if(!SymIsTerminal(config.symmetry) || config.count != params.numCatalysts) {

      LifeState activePart =
          (~history).ZOI() & config.state & ~config.startingCatalysts;

      bool hasActivePart = !activePart.IsEmpty();

      if (hasActivePart) {
        for (unsigned s = 0; s < catalysts.size(); s++) {
          if (catalysts[s].hasLocus) {
            LifeState hitLocations =
                activePart.Convolve(catalysts[s].locusAvoidMask);
            masks[s] |= hitLocations;
          }
        }
      }

      if (hasActivePart) {
        TryApplyingSymmetry(config, history, required, masks,
                            shiftedTargets, triedOffsets, missingTime,
                            recoveredTime, activePart);
      }

      // Try adding a catalyst
      if (hasActivePart && config.state.gen >= params.startGen &&
          config.count != params.numCatalysts) {
        LifeState next = config.state;
        next.Step();
        LifeState twonext = next;
        twonext.Step();

        TryAddingCatalyst(config, history, required, masks,
                          shiftedTargets, triedOffsets, missingTime,
                          recoveredTime, activePart, twonext);
        config.state = next;
        // The above also steps config.state
      } else {
        // No need to update history, not used for anything once all
        // catalysts are placed
        if (!SymIsTerminal(config.symmetry) || config.count != params.numCatalysts)
          history |= config.state;
        config.state.Step();
      }
      } else {
        config.state.Step();
      }

      for (unsigned i = 0; i < config.count; i++) {
        if (config.state.Contains(shiftedTargets[i]) || catalysts[config.curs[i]].sacrificial) {
          missingTime[i] = 0;
          recoveredTime[i] += 1;
        } else {
          missingTime[i] += 1;
          recoveredTime[i] = 0;
        }

        if (missingTime[i] > catalysts[config.curs[i]].maxDisappear) {
          failuretime = config.state.gen;
          failure = true;
          break;
        }
      }

      if (config.count == params.numCatalysts && !success) {
        bool allRecovered = true;
        for (unsigned i = 0; i < config.count; i++) {
          if (catalysts[config.curs[i]].sacrificial)
            continue;
          if (recoveredTime[i] < params.stableInterval || missingTime[i] > 0) {
            allRecovered = false;
          }
        }
        if (allRecovered) {
          success = true;
          successtime = g;
        }
      }

      if (config.count == 0 && config.symmetry == C1)
        Report();
    }

    if (!failure)
      failuretime = filterMaxGen;

    if (success)
      ReportSolution(config, successtime, failuretime);
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
