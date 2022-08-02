# CatForce
GoL Catalyst search utility based on LifeAPI library.

The main advantage of CatForce is that it doesn't make any assumptions
about the nature of the interaction. As long as the catalysts are back
in place in a given number of generations, they all could be destroyed
and reappear several times.

<!-- The torus centre is `(0, 0)` and left upper corner is `(-32, -->
<!-- -32)` and lower right corner is `(31,31)`. It has the same Y axis as -->
<!-- Golly (up is negative Y). -->

Usage
--
Run `make`, and then `./CatForce inputfile.in`. Currently only tested with `clang`.

Input File Format
--
See `examples/p83.in` etc. Some useful lists of catalysts are given in `catlists/`.

Parameters are separated by `" "` - i.e. space. Parentheses below
denote optional parameters.

| Line                  | Parameter              | Description                                                             |
|-----------------------|------------------------|-------------------------------------------------------------------------|
| `max-gen`             | `n`                    | Overall maximum generation                                              |
| `start-gen`           | `n`                    | The min gen the _first_ encounter occurs                                |
| `last-gen`            | `n`                    | The max gen the _first_ encounter occurs                                |
| `num-catalyst`        | `n`                    | The number of catalyst to place                                         |
| `num-transparent`     | `n`                    | The number of catalysts marked `transparent` that may appear            |
| `stable-interval`     | `n`                    | Gens the catalysts must remain untouched to be considered stable        |
| `search-area`         | `x y w h`              | Search area                                                             |
| `pat`                 | `rle`                  | The active pattern                                                      |
|                       | `(dx dy)`              | Offset applied to the active pattern                                    |
| `cat `                | `rle`                  | A catalyst                                                              |
|                       | `max-active`           | Number of generations in a row the catalyst may be missing              |
|                       | `dx dy`                | Offset applied centre the catalyst, typically negative                  |
|                       | `symmetries-char`      | Character denoting the symmetry of the catalyst (see below)             |
|                       | `(forbidden rle x y)`  | Forbidden pattern around the catalyst (see below)                       |
|                       | `(required rle x y)`   | Cells of a catalyst that must stay ON in every generation               |
|                       | `(mustinclude)`        | Solutions must use at least one `mustinclude` catalyst                  |
|                       | `(transparent)`        | Marked as transparent for the purposes of `num-transparent`             |
| `output`              | `filename`             | Output filename                                                         |
| `(and/or)filter`      | `genOrRange rle dx dy` | Filter that must be matched for the solution to be accepted (see below) |
| `full-report`         | `filename`             | Output filename for solutions ignoring all pattern filters              |
| `max-category-size`   | `n`                    | Maximum output row length before more solutions are dropped             |
| `fit-in-width-height` | `w h`                  | Only allow solutions where all catalysts fit in a `w` by `h` rectangle  |
| `also-required`       | `rle x y`              | Require `rle` to be present in every generation                         |
| `symmetry`            | `symmetrycode`         | Global symmetry of the entire pattern (see below)                       |

**Catalyst Symmetry**: A character specifying what transformations are
applied to the catalyst:
- `|`: mirror by y.
- `-`: mirror by x.
- `+` or `@`: all 90 degree rotations, for D2- or D2/ invariant catalysts (eg boat).
- `x`: for 180 degree symmetrical catalysts.
- `/`: diagonal mirror.
- `*`: all 8 transformations.

**Filters**: Every potential solution will be checked against the
provided filters: only solutions that match the filter will be
kept. There are two variants: `andfilter` and `orfilter`. Every
`andfilter ...` must match and at least one `orfilter ...` must match
for the solution to succeed.

**Forbidden**: Checks for `rle` in the same location as the
catalyst. If `rle` is matched in any generation, the solution is
dropped. You can see `forbidden` as a "catalyst attached filter",
i.e. the filter location is based on the catalyst location. Intended
to exclude eaters/boat-bits and similar unwanted garbage. You may have
several forbidden patterns per catalyst. See the files in `catlists`
for examples.

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


`symmetry s` automatically applies the symmetry `s` to the active
pattern and catalysts. The options for `s` are:
- D2 symmetries: `D2|` (reflect across x = 0), `D2-` (y = 0), `D2\` (y=x) and `D2/` (y=-x)
- C2 symmetries: `C2_1` (bounding box odd by odd), `C2_2` (odd by even), `C2_4` (even by even)
- C4 symmetries: `C4_1`, (odd by odd) `C4_4` (even by even)
- D4 symmetries: `D4_+1`, `D4_+2` (odd by even),`D4_+4`, `D4_x1`, and `D4_x4`
- D8 symmetries: `D8_1` and `D8_4`

The options D2 follow
[LLS](https://gitlab.com/OscarCunningham/logic-life-search), while the
rest of the symmetries follow
[Catagolue](https://catagolue.hatsya.com/census). [LifeWiki](https://conwaylife.com/wiki/Static_symmetry)
has a more in-depth explanation, though the notation there is slightly
different. (They choose `C2_2` to be even by odd not odd by even).

There are LLS-inspired alternatives for the other groups:
- D2 symmetries: as above
- C2 symmetries: `C2 C2even C2|even C2-even`
- C4 symmetries: `C4 C4even`
- D4 symmetries: `D4+ D4+|even D4+-even D4+|even D4x D4xeven`
- D8 symmetries: `D8 D8even`

The default is "odd" transformations, i.e. those that fix the cell at
(0,0); `verticaleven` and `horizontaleven` (referring to bounding box
dimensions) are recognized as equivalent to `-even` and `|even`,
respectively.
