#include "LifeAPI.h"
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <iostream>

enum CatEvalMode {
  RANDOM_SOUPS,
  FILE_SOUPS,
  ARG_SOUP,
  DIGESTS,
};
const CatEvalMode mode = FILE_SOUPS;

enum CatEvalOrder {
  SMALLEST,
  FIRST,
};
const CatEvalOrder order = SMALLEST;

const unsigned startTime = 10;
// const unsigned endTime = 100;
const unsigned endTime = 50;

// Even, for p2
// const unsigned alignment = 16;
const unsigned alignment = 8;
const unsigned fastcheck = 40;
const unsigned fastcheckinterval = 2;
const unsigned stabletime = 8;

class CatalystData {
public:
  LifeState state;
  std::vector<SymmetryTransform> transforms;
  std::vector<LifeState> orientations;
  std::vector<LifeState> orientationsMatch;
  std::vector<LifeState> reactionMasks;
  std::vector<LifeState> reactionMasks1;
  unsigned order;
  uint64_t digest;
  CatalystData(std::string rle, unsigned order);
};

std::pair<std::vector<SymmetryTransform> ,std::vector<LifeState>> Orientations(LifeState &state) {
  std::vector<SymmetryTransform> transfs;
  std::vector<LifeState> result;

  std::array<SymmetryTransform, 8> allTransforms = {
    Identity,
    ReflectAcrossXEven,
    ReflectAcrossYEven,
    Rotate90Even,
    Rotate270Even,
    Rotate180EvenBoth,
    ReflectAcrossYeqX,
    ReflectAcrossYeqNegX,
  };

  LifeState current = state;

  for(int i = 0; i < 2; i++){
  for(auto tr : allTransforms) {
    LifeState transformed = current;
    transformed.Transform(tr);


    LifeState shifted = transformed;
    std::array<int, 4> bounds = shifted.XYBounds();
    shifted.Move(-bounds[0], -bounds[1]);

    if(std::find(result.begin(), result.end(), shifted) == result.end()) {
      transfs.push_back(tr);
      result.push_back(transformed);
    }
  }
  current.Step();
  if(current == state)
    break;
  }

  return {transfs, result};
}

CatalystData::CatalystData(std::string rle, unsigned o) {
  state = LifeState::Parse(rle.c_str());
  digest = state.GetHash();
  order = o;
  std::tie(transforms, orientations) = Orientations(state);
  for(auto &o : orientations) {
    LifeState zoi = o.ZOI();
    //    orientationsZOI.push_back(zoi);

    LifeState mask = o.BigZOI();
    mask.Transform(Rotate180OddBoth);
    mask.RecalculateMinMax();
    reactionMasks.push_back(mask);


    LifeState stepped = o;
    stepped.Step();

    LifeState mask1 = stepped.BigZOI();
    mask1.Transform(Rotate180OddBoth);
    mask1.RecalculateMinMax();
    reactionMasks1.push_back(mask1);

    LifeState corona = (o | stepped).ZOI() & ~(o | stepped);
    orientationsMatch.push_back(corona | (o & stepped));
  }
}

class Perturbation {
public:
  LifeState catalyst;
  LifeState initial;
  uint64_t order;
  uint64_t souporder;
  int postPerturbed;
  uint64_t postdigest;

  LifeState contact;
  unsigned duration;
};

std::vector<Perturbation> Perturbations(CatalystData cat, LifeState pat, int souporder) {
  std::vector<Perturbation> result;
  LifeState corona = LifeState::Parse("3o$3o$3o!");
  corona.Move(-1, -1);

  LifeState margin = ~LifeState::SolidRect(-31, -31, 62, 62);

  for(int i = 0; i < cat.orientations.size(); i++) {
    LifeState &o = cat.orientations[i];
    LifeState &r = cat.reactionMasks[i];
    LifeState &r1 = cat.reactionMasks1[i];

    LifeState mask;
    LifeState orig = pat;
    LifeState current = pat;
    current.gen = 0;

    // LifeState afterfast = pat;
    // afterfast.Step(fastcheck);

    for(int g = 0; g < endTime; g++) {
      LifeState newPlacements(false);
      if(g%2 == 0)
        newPlacements = current.Convolve(r) & ~mask;
      else
        newPlacements = current.Convolve(r1) & ~mask;

      if (g < startTime) {
        mask |= newPlacements;
        current.Step();
        continue;
      }
      LifeState next = current;
      next.Step();

      while (!newPlacements.IsEmpty()) {
        auto newPlacement = newPlacements.FirstOn();
        newPlacements.Erase(newPlacement.first, newPlacement.second);

        LifeState zeropositioned = o;
        zeropositioned.Move(newPlacement.first, newPlacement.second);
        LifeState onepositioned = zeropositioned;
        onepositioned.Step();

        LifeState positioned;
        LifeState positionednext;
        if(g % 2 == 0) {
          positioned = zeropositioned;
          positionednext = onepositioned;
        } else {
          positioned = onepositioned;
          positionednext = zeropositioned;
        }

        if(!(positioned & margin).IsEmpty()) {
          mask.Set(newPlacement.first, newPlacement.second);
          continue;
        }

        LifeState positionedmatch = cat.orientationsMatch[i];
        positionedmatch.Move(newPlacement.first, newPlacement.second);

        LifeState currentwcat = current | positioned;
        currentwcat.gen = current.gen;

        LifeState wcatnext = currentwcat;
        wcatnext.Step();

        LifeState difference = (next | positionednext) ^ wcatnext;
        bool interacted = !difference.IsEmpty();
        if (!interacted)
          continue;

        LifeState contact = difference & wcatnext;

        mask.Set(newPlacement.first, newPlacement.second);

        currentwcat = current | positioned;
        currentwcat.gen = current.gen;

        bool recovered = false;
        int c;
        for(c = fastcheckinterval; c <= fastcheck; c += fastcheckinterval) {
          currentwcat.Step(fastcheckinterval);
          if ((positionedmatch & (currentwcat ^ positioned)).IsEmpty()) {
            recovered = true;
            break;
          }
        }
        if(!recovered)
          continue;

        int roughrecovery = currentwcat.gen;

        // Now find when the catalysis actually ends
        currentwcat = current | positioned;
        currentwcat.gen = current.gen;

        int j;
        for (j = 1; j < roughrecovery; j++) {
          currentwcat.Step();

          LifeState catalystpart = currentwcat.ComponentContaining(positioned, corona);

          catalystpart.Step(roughrecovery - j);
          if(j >= 2 && catalystpart == zeropositioned)
            break;
        }

        {
          // Make sure it stays recovered
          LifeState stablelookahead = currentwcat;
          stablelookahead.Step(stabletime);

          LifeState catalystpart = stablelookahead.ComponentContaining(positioned, corona);
          catalystpart.Step(roughrecovery - j - stabletime);
          if(catalystpart != zeropositioned)
            continue;
        }

        int t = g + j;
        int gap = alignment - (t % alignment);
        if(gap == alignment) gap = 0;
        t += gap;

        currentwcat.Step(t - g - j);

        LifeState catalystpart = currentwcat.ComponentContaining(zeropositioned, corona);

        if(t > roughrecovery && catalystpart != zeropositioned && catalystpart != onepositioned)
          continue;

        auto digest = (currentwcat & ~catalystpart).GetHash();

        currentwcat = current | positioned;
        currentwcat.gen = current.gen;
        LifeState active;

        LifeState currentpositioned = positioned;

        for (int j = 0; g+j < roughrecovery; j++) {
          currentwcat.Step();
          currentpositioned.Step();
          active |= (currentwcat ^ currentpositioned) & currentpositioned.ZOI();
        }

        auto cs = positioned.Components();
        if (cs.size() > 1) {
          // Check all components were used
          bool missed = false;
          bool hasTransparent = false;
          for (auto &c : cs) {
            if((c.ZOI() & active).IsEmpty()) {
              missed = true;
              break;
            }
            if((c & ~active).IsEmpty()) {
              hasTransparent = true;
            }
          }
          if(missed)
            continue;

          // Check that no components recover on their own

          bool worksAlone = false;

          for (auto &c : cs) {
            bool isHit = false;
            currentwcat = current | c;
            int presenttime = 0;
            for (int j = 0; g+j < roughrecovery; j++) {
              currentwcat.Step();
              bool present = ((currentwcat ^ c) & positionedmatch).IsEmpty();
              isHit = isHit || !present;
              if(isHit && present)
                presenttime += 1;
              if(isHit && !present)
                presenttime = 0;
              if(presenttime > stabletime) {
                worksAlone = true;
                break;
              }
            }
          }
          if(worksAlone)
            continue;
        }

        LifeState signature = contact | active;

        signature.Move(-newPlacement.first, -newPlacement.second);
        signature.Transform(TransformInverse(cat.transforms[i]));


        // if(!hasTransparent)
        //   continue;

        // if(cs.size() > 1) {
        //   for(auto &c : cs)
        //     std::cout << "  " << c.RLE() << std::endl;

        //   std::cout << (current | positioned).RLE() << std::endl;
        //   std::cout << positioned.RLE() << std::endl;
        //   std::cout << active.RLE() << std::endl;
        // }

        // bool completelyRecovered = (positionedzoi & (currentwcat ^ positioned)).IsEmpty();
        // if (!completelyRecovered)
        //   continue;

        // LifeState catalystpart = currentwcat.ComponentContaining(currentwcat & positioned, corona);
        // currentwcat &= ~catalystpart;

        // int t = g + j + sparkdie;
        // int gap = alignment - (t % alignment);
        // if(gap == alignment) gap = 0;
        // t += gap;

        // currentwcat.Step(t - g - j);

        // auto digest = currentwcat.GetHash();

        Perturbation p;
        p.catalyst = cat.state;
        p.initial = orig | zeropositioned;
        p.postPerturbed = t;
        p.postdigest = digest;
        p.order = cat.order;
        p.souporder = souporder;
        p.contact = signature;
        p.duration = c;
        result.push_back(p);
      }
      current = next;
    }
  }
  return result;
}

void Merge(std::map<uint64_t, Perturbation> &ps, std::vector<Perturbation> &news) {
  for (auto &p : news) {
    auto search = ps.find(p.postdigest);
    if (search == ps.end()) {
      ps[p.postdigest] = p;
    } else {
      switch (order) {
      case SMALLEST:
        if (search->second.catalyst.GetPop() > p.catalyst.GetPop()) {
          ps[p.postdigest] = p;
        }
        break;
      case FIRST:
        if (search->second.order > p.order) {
          ps[p.postdigest] = p;
        }
        break;
      }
    }
  }
}

void Report(int total, std::vector<CatalystData> &catalysts, std::map<uint64_t, Perturbation> &ps) {
  std::map<std::string, std::vector<Perturbation>> inverted;
  for (auto &c : catalysts) {
    LifeState key = c.state;
    key.Move(-32, -32);
    inverted[key.RLE()] = {};
  }
  for (auto &p : ps) {
    LifeState key = p.second.catalyst;
    key.Move(-32, -32);

    inverted[key.RLE()].push_back(p.second);
  }

  std::vector<std::pair<std::string, std::vector<Perturbation>>> pairs;
  for (auto itr = inverted.begin(); itr != inverted.end(); ++itr)
    pairs.push_back(*itr);

  std::sort(pairs.begin(), pairs.end(),
            [](std::pair<std::string, std::vector<Perturbation>> a,
               std::pair<std::string, std::vector<Perturbation>> b)
            { return a.second.size() > b.second.size(); });

  for (auto &p : pairs) {
    if(p.second.size() == 0) break;
    std::cout << p.first << std::endl;
  }

  for (auto &p : pairs) {
    if(p.second.size() == 0) break;

    std::cout << p.first << ": " << (float)p.second.size() / (float)total << std::endl;

    std::map<std::string, std::vector<Perturbation>> inverted;
    for (auto &p : p.second) {
      auto keypat = p.contact;
      keypat.Move(-32, -32);
      auto key = keypat.RLE() + std::to_string(p.duration);
      if (!inverted.contains(key))
        inverted[key] = {};
      inverted[key].push_back(p);
    }
    std::vector<std::pair<std::string, std::vector<Perturbation>>> pairs;
    for (auto itr = inverted.begin(); itr != inverted.end(); ++itr)
      pairs.push_back(*itr);

    std::sort(pairs.begin(), pairs.end(),
              [](std::pair<std::string, std::vector<Perturbation>> a,
                 std::pair<std::string, std::vector<Perturbation>> b)
              { return a.second.size() > b.second.size(); });

    for (auto &r : pairs) {
      std::sort(r.second.begin(), r.second.end(),
                [](Perturbation a, Perturbation b) {
                  return a.souporder < b.souporder;
                });

      LifeState toPrint = r.second[0].initial;
      std::array<int, 4> bounds = toPrint.XYBounds();
      toPrint.Move(-bounds[0], -bounds[1]);
      toPrint.Move(-32, -32);

      std::cout << "  " << toPrint.RLE() << std::endl;
    }
  }
}

void ReportArg(int total, std::vector<CatalystData> &catalysts, std::map<uint64_t, Perturbation> &ps) {
  std::map<std::string, std::vector<LifeState>> inverted;
  for (auto &c : catalysts) {
    LifeState key = c.state;
    key.Move(-32, -32);
    inverted[key.RLE()] = {};
  }
  for (auto &p : ps) {
    LifeState key = p.second.catalyst;
    key.Move(-32, -32);

    LifeState value = p.second.initial;
    //value.Step(p.second.postPerturbed);
    std::array<int, 4> bounds = value.XYBounds();
    value.Move(-bounds[0], -bounds[1]);
    value.Move(-32, -32);

    inverted[key.RLE()].push_back(value);
  }

  std::vector<std::pair<std::string, std::vector<LifeState>>> pairs;
  for (auto itr = inverted.begin(); itr != inverted.end(); ++itr)
    pairs.push_back(*itr);

  std::sort(pairs.begin(), pairs.end(),
            [](std::pair<std::string, std::vector<LifeState>> a,
               std::pair<std::string, std::vector<LifeState>> b)
            { return a.second.size() > b.second.size(); });

  for (auto &p : pairs) {
    for (auto &r : p.second) {
      std::cout << r.RLE() << std::endl;
    }
  }
}

int main(int argc, char *argv[]) {
  char *fname = argv[1];
  std::ifstream infile;
  infile.open(fname, std::ifstream::in);

  std::vector<CatalystData> catalysts;

  std::string line;
  unsigned order = 0;
  while (std::getline(infile, line)) {
    catalysts.push_back(CatalystData(line, order));
    order++;
  }

  std::vector<std::string> filesoups;
  if (mode == FILE_SOUPS) {
    char *soupfname = argv[2];
    std::ifstream soupfile;
    soupfile.open(soupfname, std::ifstream::in);

    std::string line;
    while (std::getline(soupfile, line)) {
      filesoups.push_back(line);
    }
  }

  std::map<uint64_t, Perturbation> ps;
  int end = 10000;
  if(mode == FILE_SOUPS)
    end = filesoups.size();

  for (int i = 0; i < end; i++) {
    LifeState soup;
    switch(mode) {
    case RANDOM_SOUPS: {
      soup = LifeState::RandomState() & LifeState::SolidRect(0, 0, 12, 12);
      std::cout << "Doing soup " << soup.RLE() << std::endl;
      break;
    }
    case FILE_SOUPS: {
      soup = LifeState::Parse(filesoups[i].c_str());
      std::cout << "Doing soup " << soup.RLE() << std::endl;
      break;
    }
    case ARG_SOUP:
      soup = LifeState::Parse(argv[2]);
      break;
    case DIGESTS: {
      soup = LifeState::RandomState() & LifeState::SolidRect(0, 0, 12, 12);
      break;
    }
    }

    // soup = LifeState::Parse("2o3b4obo$b3ob2ob4o$2b2obo5bo$5b3o$b2obo2b2obo$4b3o2bo$bo2bobobob2o$bo3bob5o$2bo2b2o2bo$ob3o2bobo$4ob2obo2bo$o2b3obo!");

    uint64_t souphash = soup.GetHash();
    for (auto &c : catalysts) {
      auto news = Perturbations(c, soup, i);

      if(mode == DIGESTS) {
        if(news.size() > 0) {
          std::cout << c.state.RLE() << " ";

          std::vector<uint64_t> seenhashes;
          for(auto &p : news) {
            uint64_t longhash = (p.postdigest ^ souphash);
            if(std::find(seenhashes.begin(), seenhashes.end(), longhash) == seenhashes.end()) {
              seenhashes.push_back(longhash);
              std::cout << longhash << " ";
            }
          }
          std::cout << std::endl;
        }
      }
      else
        Merge(ps, news);
    }

    // Report(i, catalysts, ps);
    // exit(0);

    switch(mode) {
    case RANDOM_SOUPS:
    case FILE_SOUPS: {
      if(i%10 == 0)
        Report(i+1, catalysts, ps);
      break;
    }
    case ARG_SOUP: {
      ReportArg(i+1, catalysts, ps);
      exit(0);
    }
    }

  }
  Report(end, catalysts, ps);
  exit(0);
}
