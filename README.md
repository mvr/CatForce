# CatForce
GOL Catalyst search utility based on LifeAPI library.

The main advantage of CatForce is that it doesn't make any assumptions
about the nature of the interaction. As long as the catalysts are back
in place in a given number of generations, they all could be destroyed
and reappear several times.

It uses brute force search in an area rather than a tree search.

<!-- The torus centre is `(0, 0)` and left upper corner is `(-32, -->
<!-- -32)` and lower right corner is `(31,31)`. It has the same Y axis as -->
<!-- Golly (up is negative Y). -->

<!-- Another feature of LifeAPI preserved in CatForce is eliminating edge -->
<!-- gliders. LifeAPI is currently "listening" on the edges of the torus -->
<!-- for gliders and removes them. This feature is also true for CatForce. -->

Usage
--
Run `make`, and then `./CatForce inputfile.in`. Currently only tested with `clang`.

Input File Format
--
See `1.in` for an example.

The option delimiter is `" "` - i.e. space. 

| Line                        | Parameter                | Description                                                            |
|-----------------------------|--------------------------|------------------------------------------------------------------------|
| `max-gen`                   | `n`                      | Overall maximum generation                                             |
| `start-gen`                 | `n`                      | The min gen the first encounter must happen by                         |
| `last-gen`                  | `n`                      | The max gen the _first_ encounter is allowed                           |
| `num-catalyst`              | `n`                      | The number of catalyst to place                                        |
| `stable-interval`           | `n`                      | Gens the catalysts must remain untouched to be considered stable       |
| `search-area`               | `x y w h`                | Search area                                                            |
| `pat`                       | `rle`                    | The active pattern                                                     |
|                             | `(dx dy)`                | Optional offset                                                        |
| `cat `                      | `rle`                    | A catalyst                                                             |
|                             | `max-active`             | Number of generations in a row the catalyst may be missing             |
|                             | `dx dy`                  | Offset to centre the catalyst (typically negative)                     |
|                             | `symmetries-char`        | Character denoting the symmetry of the catalyst                        |
|                             | `(forbidden rle x y)`    | Optional forbidden pattern around the catalyst                         |
| `output`                    | `filename`               | Output filename                                                        |
| `filter`                    | `gen rle dx dy`          | Filter that must be matched for the solution to be accepted            |
| `filter`                    | `min-max rle dx dy`      | Filter in a range of generations.                                      |
| `filter`                    | `min-max rle dx dy sym`  | Filter with symmetry (see below)                                       |
| `full-report`               | `filename`               | Output filename for solutions ignoring all pattern filters             |
| `max-category-size`         | `n`                      | Maximum output row length before more solutions are dropped            |
| `fit-in-width-height`       | `w h`                    | Only allow solutions where all catalysts fit in a `w` by `h` rectangle |
| `symmetry`                  |                          | See below                                                              |
| `stop-after-cats-destroyed` | `n`                      | Filters must be met n generations after catalyst destruction or sooner |
| `also-required`             | `rle dx dy`              | Require pattern to be present in every generation.                     |
| `quiet-mode`                | `true`                   | Omit "placed catalyst [] at []" messages.                              |
| `report-matches`            | `true`                   | At end of search print what generation matched any gen range filters.  |
| `first-transparent`         | `true`                   | Only search reactions that have early-on transparent catalyses.        |



**Encounter**: An active cell from the input pattern is present in
the immediate neighbourhood of the catalyst.

**Catalyst Symmetry**: A character specifying what transformations are
applied to the catalyst:
- `|` mirror by y.
- `-` mirror by x.
- `+` all 90 degree rotations, for D2- or D2/ invariant catalysts (eg boat).
- `/` diagonal mirror
- `x` for 180 degree symmetrical catalysts.
- `*` all 8 transformations.

**Forbidden**: Will search for `rle` in the same location as the
catalyst. If `rle` is matched in any generation, the solution is
omitted from the report. You can see forbidden as "catalyst attached
filter", i.e. the forbidden pattern is based on catalyst location.
Made to exclude eaters/boat-bits and similar unwanted "garbage". You
may have several forbidden patterns per catalyst. See `3.in` for
an example.

**Filters**: one can use several filters. Every filter will be
checked if successful catalyst was found. Each filter will select
for patterns matching the live cells in `rle` and dead cells in close
proximity neighbourhood to the live cells in `rle`. Filters also
accept a symmetry group: one of the transformed images of the
`rle` must appear. See the symmetric searches section for syntax.
(Intended for searching for oscillators that flip every half period:
 this way, one search and one filter covers both re-occurrance
 and the various the oscillator could flip.)

**Stop After Catalysts Destroyed**: if the catalysts all recover,
only check filters out to n generations after the first
catalyst is destroyed. (This cuts down on false positives 
when searching for oscillators.)

<!-- Combining Results -->
<!-- --- -->

<!-- `combine-results yes [<survive-0> <survive-1> ...]` -->

<!-- If this feature is enabled the search will at first ignore all filters -->
<!-- and survival inputs, and will search all the possible catalysts. Then -->
<!-- it will try to combine all the found catalysts in all possible -->
<!-- combinations, and only then will filter by `survive-i` and apply the -->
<!-- filters to exclude them from the final report. -->
 
<!-- This feature will generate report as follows: -->
 
<!-- - `output.rle` - all the possible catalysts. -->
<!-- - `output.rle_Combined*.rle` - will generate all combined reports. -->
<!-- - `output.rle_Final.rle` - the final report. **This is the main output.** -->
 
<!-- Optional survival filter per "iteration" are added. Combine works as -->
<!-- follows: each time it start from the initial search results (combine -->
<!-- by default uses survive count = 1), and tries to add catalyst from -->
<!-- those results. Sometimes one could get explosion, if the interaction -->
<!-- is very potent. So filter is added to limit the combine, by surviving -->
<!-- count (if something doesn't survive with two catalyst for 5 -->
<!-- iterations, it's probably junk - so CatForce will filter it on the -->
<!-- second combine iteration and not in the end). -->
 
<!-- This allows faster and more efficient combine operation with very -->
<!-- potent conduits which otherwise would overflow the system, with many -->
<!-- useless catalysts. -->
 
<!-- **NOTE** Recommended for use only for `num-catalyst` = 1 or 2 -->

<!-- **NOTE** See 4.in file for example.  -->
 
<!-- **NOTE** CatForce will use the last `survive-i` as the default from -->
<!-- that point on. If you don't enter any numbers it will use survival -->
<!-- count 1, and will filter only when finish all possible combinations. -->

Symmetric Searches
---

`symmetry s`

Automatically apply the symmetry `s` to the active pattern and catalysts. Options:
- D2 symmetries: `D2|` (reflect across x = 0), `D2-` (y = 0), `D2\` (y=x) and `D2/` (y=-x)
- C2 symmetries: `C2_1` (bounding box odd by odd), `C2_2` (odd by even), `C2_4` (even by even)
- C4 symmetries: `C4_1`, (odd by odd) `C4_4` (even by even)
- D4 symmetries: `D4_+1`, `D4_+2` (odd by even),`D4_+4`, `D4_x1`, and `D4_x4`
- D8 symmetries: `D8_1` and `D8_4`


The options D2 follow [LLS](https://gitlab.com/OscarCunningham/logic-life-search), while the rest of the symmetries follow [Catagolue](https://catagolue.hatsya.com/census). [LifeWiki](https://conwaylife.com/wiki/Static_symmetry) has a more in-depth explanation, though the notation there is slightly
different. (They choose `C2_2` to be even by odd not odd by even).

There's LLS-inspired alternatives for the other groups: 
- D2 symmetries: as above
- C2 symmetries: `C2 C2even C2|even C2-even`
- C4 symmetries: `C4 C4even`
- D4 symmetries: `D4+ D4+|even D4+-even D4+|even D4x D4xeven`
- D8 symmetries: `D8 D8even`

Here, things default to odd transformations, ie ones that fix the cell at (0,0);`verticaleven` and `horizontaleven` (referring to bounding
box dimensions) are recognized as equivalent to `-even` and `|even`,
respectively.

