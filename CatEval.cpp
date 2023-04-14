#include "LifeAPI.h"
#include <vector>
#include <map>
#include <fstream>
#include <iostream>

// Even, for p2
const unsigned alignment = 16;
const unsigned fastcheck = 30;
const unsigned fastcheckinterval = 6;
const unsigned stabletime = 8;

class CatalystData {
public:
  LifeState state;
  std::vector<LifeState> orientations;
  std::vector<LifeState> orientationsZOI;
  std::vector<LifeState> orientationsMatch;
  std::vector<LifeState> reactionMasks;
  std::vector<LifeState> reactionMasks1;
  CatalystData(std::string rle);
};

std::vector<LifeState> Orientations(LifeState &state) {
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

    std::array<int, 4> bounds = transformed.XYBounds();
    transformed.Move(-bounds[0], -bounds[1]);

    if(std::find(result.begin(), result.end(), transformed) == result.end()) {
      result.push_back(transformed);
    }
  }
  current.Step();
  if(current == state)
    break;
  }

  return result;
}

CatalystData::CatalystData(std::string rle) {
  state = LifeState::Parse(rle.c_str());
  orientations = Orientations(state);
  for(auto &o : orientations) {
    LifeState zoi = o.ZOI();
    orientationsZOI.push_back(zoi);

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
  int postPerturbed;
  uint64_t postdigest;
};

std::vector<Perturbation> Perturbations(CatalystData cat, LifeState pat) {
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

    for(int g = 0; g < 50; g++) {
      LifeState newPlacements(false);
      if(g%2 == 0)
        newPlacements = current.Convolve(r) & ~mask;
      else
        newPlacements = current.Convolve(r1) & ~mask;

      if (g < 10) {
        mask |= newPlacements;
        current.Step();
        continue;
      }
      LifeState next = current;
      next.Step();
      LifeState next2 = next;
      next2.Step();

      while (!newPlacements.IsEmpty()) {
        auto newPlacement = newPlacements.FirstOn();
        newPlacements.Erase(newPlacement.first, newPlacement.second);

        LifeState zeropositioned = o;
        zeropositioned.Move(newPlacement.first, newPlacement.second);
        LifeState onepositioned = zeropositioned;
        onepositioned.Step();

        LifeState positioned;
        if(g % 2 == 0)
          positioned = zeropositioned;
        else
          positioned = onepositioned;

        if(!(positioned & margin).IsEmpty()) {
          mask.Set(newPlacement.first, newPlacement.second);
          continue;
        }

        // LifeState positionedzoi = z;
        // positionedzoi.Move(newPlacement.first, newPlacement.second);
        // LifeState positionedbigzoi = bz;
        // positionedbigzoi.Move(newPlacement.first, newPlacement.second);
        // LifeState positionedcorona = c;
        // positionedcorona.Move(newPlacement.first, newPlacement.second);

        LifeState positionedmatch = cat.orientationsMatch[i];
        positionedmatch.Move(newPlacement.first, newPlacement.second);

        LifeState currentwcat = current | positioned;
        currentwcat.gen = current.gen;

        LifeState wcatnext = currentwcat;
        wcatnext.Step();
        wcatnext.Step();

        bool interacted = !((next2 | positioned) ^ wcatnext).IsEmpty();
        if (!interacted)
          continue;

        mask.Set(newPlacement.first, newPlacement.second);

        currentwcat = current | positioned;
        currentwcat.gen = current.gen;

        bool recovered = false;
        int c;
        for(c = 0; c <= fastcheck; c += fastcheckinterval) {
          currentwcat.Step(fastcheckinterval);
          if ((positionedmatch & (currentwcat ^ positioned)).IsEmpty()) {
            recovered = true;
            break;
          }
        }
        if(!recovered)
          continue;

        int roughrecovery = currentwcat.gen;

        // Make sure it stays recovered
        currentwcat.Step(stabletime);
        if (!(positionedmatch & (currentwcat ^ positioned)).IsEmpty()) {
          continue;
        }

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

        int t = g + j;
        int gap = alignment - (t % alignment);
        if(gap == alignment) gap = 0;
        t += gap;

        currentwcat.Step(t - g - j);


        LifeState catalystpart = currentwcat.ComponentContaining(zeropositioned, corona);

        if(t > roughrecovery+stabletime && catalystpart != zeropositioned && catalystpart != onepositioned)
          continue;

        auto digest = (currentwcat & ~catalystpart).GetHash();

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
      if (search->second.catalyst.GetPop() > p.catalyst.GetPop()) {
        ps[p.postdigest] = p;
      }
    }
  }
}

bool CompareLength(std::pair<std::string, std::vector<LifeState>> &a, std::pair<std::string, std::vector<LifeState>> &b) { return a.second.size() > b.second.size();}

void Report(int total, std::vector<CatalystData> &catalysts, std::map<uint64_t, Perturbation> &ps) {
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

  sort(pairs.begin(), pairs.end(), CompareLength);

  for (auto &p : pairs) {
    //    if(p.second.size() == 0) break;
    std::cout << p.first << std::endl;
  }

  for (auto &p : pairs) {
    //    if(p.second.size() == 0) break;
    std::cout << p.first << ": " << (float)p.second.size() / (float) total << std::endl;
    int cnt = 0;
    for (auto &r : p.second) {
      std::cout << "  " << r.RLE() << std::endl;
      //      if (cnt > 50)
      if (cnt > 10)
        break;
      cnt++;
    }
  }
}

int main(int argc, char *argv[]) {
  char *fname = argv[1];
  std::ifstream infile;
  infile.open(fname, std::ifstream::in);

  std::vector<CatalystData> catalysts;

  std::string line;
  while (std::getline(infile, line)) {
    catalysts.push_back(CatalystData(line));
  }

  std::map<uint64_t, Perturbation> ps;

  for (int i = 1; i <= 10000; i++) {
    LifeState soup =
        LifeState::RandomState() & LifeState::SolidRect(0, 0, 12, 12);
    std::cout << "Doing soup " << soup.RLE() << std::endl;

    // LifeState soup = LifeState::Parse("bob4o$6o2bo2bo$4ob3o2bo$o3b4obobo$bob2o2bo2b2o$2ob2obob2o$b2ob2o2bo$bo2b4o3bo$2o3bo2bo2bo$bo3bo3bo$2bob2o3b2o$o3b5obo!");

    uint64_t souphash = soup.GetHash();
    for (auto &c : catalysts) {
      auto news = Perturbations(c, soup);

      // if(news.size() > 0) {
      //   std::cout << c.state.RLE() << " ";

      //   std::vector<uint64_t> seenhashes;
      //   for(auto &p : news) {
      //     uint64_t longhash = (p.postdigest ^ souphash);
      //     if(std::find(seenhashes.begin(), seenhashes.end(), longhash) == seenhashes.end()) {
      //       seenhashes.push_back(longhash);
      //       std::cout << longhash << " ";
      //     }
      //   }
      //   std::cout << std::endl;
      // }

      Merge(ps, news);
    }
    if(i%100 == 0)
      Report(i, catalysts, ps);
  }

  exit(0);
}
