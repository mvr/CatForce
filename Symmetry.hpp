#pragma once

#include "LifeAPI.hpp"

enum SymmetryTransform : uint32_t {
  Identity,
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
  // reflect across y = -x+3/2, fixing (0,0), instead of y=-x+1/2,
  // sending (0,0) to (-1,-1). Needed for D4x_1 symmetry.
  ReflectAcrossYeqNegXP1
};

static SymmetryTransform allTransforms[] = {
    Identity,
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
    ReflectAcrossYeqNegXP1
};

SymmetryTransform TransformInverse(SymmetryTransform transf) {
  switch (transf) {
  case Rotate90Even: return Rotate270Even;
  case Rotate90: return Rotate270;
  case Rotate270Even: return Rotate90Even;
  case Rotate270: return Rotate90;
  default: return transf;
  }
}

enum StaticSymmetry : uint32_t {
  C1,
  D2AcrossX,
  D2AcrossXEven,
  D2AcrossY,
  D2AcrossYEven,
  D2negdiagodd,
  D2diagodd,
  C2,
  C2even,
  C2verticaleven,
  C2horizontaleven,
  C4,
  C4even,
  D4,
  D4even,
  D4verticaleven,
  D4horizontaleven,
  D4diag,
  D4diageven,
  D8,
  D8even,
};

static StaticSymmetry allSymmetries[] = {
    C1,
    D2AcrossX,
    D2AcrossXEven,
    D2AcrossY,
    D2AcrossYEven,
    D2negdiagodd,
    D2diagodd,
    C2,
    C2even,
    C2verticaleven,
    C2horizontaleven,
    C4,
    C4even,
    D4,
    D4even,
    D4verticaleven,
    D4horizontaleven,
    D4diag,
    D4diageven,
    D8,
    D8even,
};



void LifeState::JoinWSymChain(const LifeState &state, int x, int y,
                   const std::vector<SymmetryTransform> &symChain) {
  // instead of passing in the symmetry group {id, g_1, g_2,...g_n} and
  // applying each to default orientation we pass in a "chain" of symmetries
  // {h_1, ...h_n-1} that give the group when "chained together": g_j =
  // product of h_1 thru h_j that way, we don't need to initialize a new
  // LifeState for each symmetry.

  LifeState transformed = state.Moved(x, y);

  for (auto sym : symChain) {
    transformed |= transformed.Transformed(sym);
  }
  *this |= transformed;
}

void LifeState::JoinWSymChain(const LifeState &state,
                   const std::vector<SymmetryTransform> &symChain) {
  LifeState transformed = state;

  for (auto sym : symChain) {
    transformed |= transformed.Transformed(sym);
  }
  *this |= transformed;
}

void LifeState::Transform(SymmetryTransform transf) {
  switch (transf) {
  case Identity:
    break;
  case ReflectAcrossXEven:
    FlipX();
    break;
  case ReflectAcrossX:
    FlipX();
    Move(0, 1);
    break;
  case ReflectAcrossYEven:
    FlipY();
    break;
  case ReflectAcrossY:
    FlipY();
    Move(1, 0);
    break;
  case Rotate180EvenBoth:
    FlipX();
    FlipY();
    break;
  case Rotate180EvenVertical:
    FlipX();
    FlipY();
    Move(1, 0);
    break;
  case Rotate180EvenHorizontal:
    FlipX();
    FlipY();
    Move(0, 1);
    break;
  case Rotate180OddBoth:
    FlipX();
    FlipY();
    Move(1, 1);
    break;
  case ReflectAcrossYeqX:
    Transpose(false);
    break;
  case ReflectAcrossYeqNegX:
    Transpose(true);
    break;
  case ReflectAcrossYeqNegXP1:
    Transpose(true);
    Move(1, 1);
    break;
  case Rotate90Even:
    FlipX();
    Transpose(false);
    break;
  case Rotate90:
    FlipX();
    Transpose(false);
    Move(1, 0);
    break;
  case Rotate270Even:
    FlipY();
    Transpose(false);
    break;
  case Rotate270:
    FlipY();
    Transpose(false);
    Move(0, 1);
    break;
  }
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

constexpr LifeState FundamentalDomain(const StaticSymmetry sym) {
  switch (sym) {
  case StaticSymmetry::C1:
    return LifeState::ConstantParse(
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$"
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$"
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$"
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o!");
  case StaticSymmetry::D2AcrossY:
  case StaticSymmetry::D2AcrossYEven:
    return LifeState::ConstantParse(
        "33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$"
        "33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$"
        "33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$"
        "33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o!");
  case StaticSymmetry::D2AcrossX:
  case StaticSymmetry::D2AcrossXEven:
    return LifeState::ConstantParse(
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$"
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o!");
  case StaticSymmetry::D2diagodd:
    return LifeState::ConstantParse(
        "2o$3o$4o$5o$6o$7o$8o$9o$10o$11o$12o$13o$14o$15o$16o$17o$18o$19o$20o$"
        "21o$22o$23o$24o$25o$26o$27o$28o$29o$30o$31o$32o$33o$34o$35o$36o$37o$"
        "38o$39o$40o$41o$42o$43o$44o$45o$46o$47o$48o$49o$50o$51o$52o$53o$54o$"
        "55o$56o$57o$58o$59o$60o$61o$62o$63o$64o$64o!");
  case StaticSymmetry::D2negdiagodd:
    return LifeState::ConstantParse(
        "64o$64o$64o$63o$62o$61o$60o$59o$58o$57o$56o$55o$54o$53o$52o$51o$50o$"
        "49o$48o$47o$46o$45o$44o$43o$42o$41o$40o$39o$38o$37o$36o$35o$34o$33o$"
        "32o$31o$30o$29o$28o$27o$26o$25o$24o$23o$22o$21o$20o$19o$18o$17o$16o$"
        "15o$14o$13o$12o$11o$10o$9o$8o$7o$6o$5o$4o$3o!");
  case StaticSymmetry::C2:
  case StaticSymmetry::C2even:
  case StaticSymmetry::C2horizontaleven:
  case StaticSymmetry::C2verticaleven:
    return LifeState::ConstantParse(
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$"
        "64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o$64o!");
  case StaticSymmetry::C4:
  case StaticSymmetry::C4even:
  case StaticSymmetry::D4:
  case StaticSymmetry::D4even:
  case StaticSymmetry::D4horizontaleven:
  case StaticSymmetry::D4verticaleven:
    return LifeState::ConstantParse(
        "33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$"
        "33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o$33o!");
  case StaticSymmetry::D4diag:
  case StaticSymmetry::D4diageven:
    return LifeState::ConstantParse(
        "2o$3o$4o$5o$6o$7o$8o$9o$10o$11o$12o$13o$14o$15o$16o$17o$18o$19o$20o$"
        "21o$22o$23o$24o$25o$26o$27o$28o$29o$30o$31o$32o$33o$34o$33o$32o$31o$"
        "30o$29o$28o$27o$26o$25o$24o$23o$22o$21o$20o$19o$18o$17o$16o$15o$14o$"
        "13o$12o$11o$10o$9o$8o$7o$6o$5o$4o$3o!");
  case StaticSymmetry::D8:
  case StaticSymmetry::D8even:
    return LifeState::ConstantParse(
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
    sym |= state;
    return sym;
  }
  case C4: {
    LifeState sym = state;
    sym.Transform(Rotate90);
    sym.Move(offset);
    sym |= state;

    LifeState sym2 = sym;
    sym2.FlipX();
    sym2.FlipY();
    sym2.Move(offset.first - offset.second + 1, offset.second + offset.first + 1);
    sym |= sym2;
    return sym;
  }
  case D2AcrossX: {
    LifeState sym = state;
    sym.FlipX();
    sym.Move(offset.first, offset.second + 1);
    sym |= state;
    return sym;
  }
  case D2AcrossY: {
    LifeState sym = state;
    sym.FlipY();
    sym.Move(offset.first + 1, offset.second);
    sym |= state;
    return sym;
  }
  case D2diagodd: {
    LifeState sym = state;
    sym.Transform(ReflectAcrossYeqX);
    sym.Move(offset);
    sym |= state;
    return sym;
  }
  case D2negdiagodd: {
    LifeState sym = state;
    sym.Transpose(true);
    sym.Move(offset.first + 1, offset.second + 1);
    sym |= state;
    return sym;
  }

  case D4: {
    LifeState acrossx = state;
    auto xoffset = PerpComponent(ReflectAcrossX, offset);
    acrossx.FlipX();
    acrossx.Move(xoffset.first, xoffset.second + 1);
    acrossx |= state;

    LifeState acrossy = acrossx;
    auto yoffset = PerpComponent(ReflectAcrossY, offset);
    acrossy.FlipY();
    acrossy.Move(yoffset.first + 1, yoffset.second);
    acrossy |= acrossx;

    return acrossy;
  }
  case D4diag: {
    LifeState acrossy = state;
    acrossy.Transform(ReflectAcrossYeqX);
    auto yoffset = PerpComponent(ReflectAcrossYeqX, offset);
    acrossy.Move(yoffset);
    acrossy |= state;

    LifeState acrossx = acrossy;
    acrossx.Transpose(true);
    auto xoffset = PerpComponent(ReflectAcrossYeqNegXP1, offset);
    acrossx.Move(xoffset.first + 1, xoffset.second + 1);
    acrossx |= acrossy;

    return acrossx;
  }

  default:
    __builtin_unreachable();
  }
}

// On intel there is a single instruction for this
// Taken from Hacker's Delight
uint64_t compress_right(uint64_t x, uint64_t m) {
   uint64_t mk, mp, mv, t;
   int i;

   x = x & m;           // Clear irrelevant bits.
   mk = ~m << 1;        // We will count 0's to right.

   for (i = 0; i < 6; i++) {
      mp = mk ^ (mk << 1);             // Parallel prefix.
      mp = mp ^ (mp << 2);
      mp = mp ^ (mp << 4);
      mp = mp ^ (mp << 8);
      mp = mp ^ (mp << 16);
      mp = mp ^ (mp << 32);
      mv = mp & m;                     // Bits to move.
      m = (m ^ mv) | (mv >> (1 << i));   // Compress m.
      t = x & mv;
      x = (x ^ t) | (t >> (1 << i));     // Compress x.
      mk = mk & ~mp;
   }
   return x;
}

inline LifeState LifeState::Halve() const {
  LifeState result;
  for(int i = 0; i < N/2; i++){
    uint64_t halvedColumn = compress_right(state[2*i], 0x5555555555555555ULL);
    halvedColumn |= halvedColumn << N/2;
    result.state[i] = halvedColumn;
    result.state[i + N/2] = halvedColumn;
  }
  return result;
}

inline LifeState LifeState::HalveX() const {
  LifeState result;
  for(int i = 0; i < N/2; i++){
    result.state[i] = state[2*i];
    result.state[i + N/2] = state[2*i];
  }
  return result;
}

inline LifeState LifeState::HalveY() const {
  LifeState result;
  for(int i = 0; i < N; i++){
    uint64_t halvedColumn = compress_right(state[i], 0x5555555555555555ULL);
    halvedColumn |= halvedColumn << N/2;
    result.state[i] = halvedColumn;
  }
  return result;
}

// (x, y) |-> (x, y + x)
inline LifeState LifeState::Skew() const {
  LifeState result;
  for(int i = 0; i < N; i++){
    result.state[i] = std::rotl(state[i], i);
  }
  return result;
}

// (x, y) |-> (x, y - x)
inline LifeState LifeState::InvSkew() const {
  LifeState result;
  for(int i = 0; i < N; i++){
    result.state[i] = std::rotr(state[i], i);
  }
  return result;
}

  LifeState IntersectingOffsets(const LifeState &pat1, const LifeState &pat2, StaticSymmetry sym) {
    LifeState transformed;
    switch (sym) {
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

  LifeState IntersectingOffsets(const LifeState &active, StaticSymmetry sym) {
    return IntersectingOffsets(active, active, sym);
  }

uint64_t LifeState::GetOctoHash() const {
    uint64_t result = 0;

    auto allTransforms = {
        Identity,           ReflectAcrossXEven,   ReflectAcrossYeqX,
        ReflectAcrossYEven, ReflectAcrossYeqNegX, Rotate90Even,
        Rotate270Even,      Rotate180EvenBoth};

    for (auto t : allTransforms) {
      LifeState transformed = Transformed(t);
      auto [x, y,_x2, y2] = transformed.XYBounds();
      transformed.Move(-x, -y);
      result ^= transformed.GetHash();
    }

    return result;
  }

LifeState LifeState::MatchesLiveAndDeadSym(const LifeState &live, const LifeState &dead) const {
    auto allTransforms = {
        Identity,           ReflectAcrossXEven,   ReflectAcrossYeqX,
        ReflectAcrossYEven, ReflectAcrossYeqNegX, Rotate90Even,
        Rotate270Even,      Rotate180EvenBoth};

    LifeState result;

    for (auto t : allTransforms) {
      LifeState matches = MatchLiveAndDead(live.Transformed(t), dead.Transformed(t));
      result |= matches.Convolve(live.Transformed(t));
    }

    return result;
  }
