# CatForce
GoL Catalyst search utility based on LifeAPI library.

<!-- The torus centre is `(0, 0)` and left upper corner is `(-32, -->
<!-- -32)` and lower right corner is `(31,31)`. It has the same Y axis as -->
<!-- Golly (up is negative Y). -->

Usage
--
Run `make`, and then `./CatForce inputfile.in`. Currently only tested with `clang`.

Input File Format
--
See `examples/p83.in` etc. Some useful lists of catalysts are given in `catlists/`.

Parameters are separated by `" "` - i.e. space. Parameters in
parentheses are optional.

| Line                  | Parameter                | Description                                                             |
|-----------------------|--------------------------|-------------------------------------------------------------------------|
| `max-gen`             | `n`                      | Overall maximum generation                                              |
| `start-gen`           | `n`                      | The min gen the _first_ encounter occurs                                |
| `last-gen`            | `n`                      | The max gen the _first_ encounter occurs                                |
| `num-catalyst`        | `n`                      | The number of catalyst to place                                         |
| `num-transparent`     | `n`                      | The number of catalysts marked `transparent` that may appear            |
| `stable-interval`     | `n`                      | Gens the catalysts must remain untouched to be considered stable        |
| `search-area`         | `x y w h`                | Search area                                                             |
| `pat`                 | `rle`                    | The active pattern                                                      |
|                       | `(dx dy)`                | Offset applied to the active pattern                                    |
| `cat `                | `rle`                    | A catalyst                                                              |
|                       | `max-active`             | Number of generations in a row the catalyst may be missing              |
|                       | `dx dy`                  | Offset applied centre the catalyst, typically negative                  |
|                       | `symmetry-char`          | Character denoting the symmetry of the catalyst (see below)             |
|                       | `(forbidden rle x y)`    | Forbidden pattern around the catalyst (see below)                       |
|                       | `(required rle x y)`     | Cells of a catalyst that must stay ON in every generation               |
|                       | `(antirequired rle x y)` | Cells of a catalyst that must stay OFF in every generation              |
|                       | `(locus rle x y)`        | Cells of a catalyst that must interact first                            |
|                       | `(must-include)`         | Solutions must use at least one `mustinclude` catalyst                  |
|                       | `(transparent)`          | Marked as transparent for the purposes of `num-transparent`             |
|                       | `(check-recovery)`       | Always check the catalyst is recovered after exactly `max-active` gens  |
|                       | `(sacrificial)`          | Does not need to recover                                                |
|                       | `(fixed)`                | Must be placed exactly at the origin                                    |
| `output`              | `filename`               | Output filename                                                         |
| `(and/or)filter`      | `genOrRange rle dx dy`   | Filter that must be matched for the solution to be accepted (see below) |
| `match`               | `genOrRange rle`         | Filter that checks whether the pattern occurs anywhere                  |
| `max-junk`            | `n`                      | For filters, the maximum non-matching non-catalyst population           |
| `full-report`         | `filename`               | Output filename for solutions ignoring all pattern filters              |
| `max-category-size`   | `n`                      | Maximum output row length before more solutions are dropped             |
| `fit-in-width-height` | `w h`                    | Only allow solutions where all catalysts fit in a `w` by `h` rectangle  |
| `also-required`       | `rle x y`                | Require `rle` to be present in every generation                         |
| `symmetry`            | `symmetry-code`          | Global symmetry of the entire pattern (see below)                       |

**Catalyst Symmetry**: A character specifying what transformations are
applied to the catalyst:
- `|`: reflect across y, for D4x catalysts (e.g. ship)
- `-`: reflect across x
- `/`: reflect diagonally, for D4+ catalysts (e.g. table on table)
- `\`: reflect diagonally the other way
- `x`: reflect diagonally both ways, for C2 symmetrical catalysts (e.g. snake)
- `+` or `@`: rotate all 90 degree increments, for D2- or D2/ catalysts (e.g. boat)
- `*`: all transforms, for C1 catalysts

**Filters**: Every potential solution will be checked against the
provided filters: only solutions that match the filters will be
kept. There are two variants: `andfilter` and `orfilter`. Every
`andfilter ...` must match and at least one `orfilter ...` must match
for the solution to succeed. A `match` filter pattern-matches the
solution for the specified pattern, and succeeds if the pattern occurs
anywhere (and otherwise behaves like an `andfilter`).

Filters that use a range of generations will only succeed after every
catalyst has interacted.

**Forbidden**: Checks for `rle` in the same location as the
catalyst. If `rle` is matched in any generation, the solution is
dropped. You can see `forbidden` as a "catalyst attached filter",
i.e. the filter location is based on the catalyst location. Intended
to exclude eaters/boat-bits and similar unwanted garbage. You may have
several forbidden patterns per catalyst. See the files in `catlists`
for examples.

**Check Recovery**: Catalysts with the `check-recovery` attribute are
immediately tested for whether they recover in `max-active`
generations, without the support of any further catalysts. This is
useful for catalysts with a bait still life that has a long recovery
time, like the hive-pushes or the loaf-spin catalysts.

Symmetric Searches
---

`symmetry s` automatically applies the symmetry `s` to the active
pattern and catalysts. The options for `s` are:
- D2 symmetries: `D2|` (reflect across x = 0), `D2|even` (x=-0.5), `D2-` (y=0), `D2-even` (y=-0.5), `D2\` (y=x) and `D2/` (y=-x)
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

Running Tests
---

You will need to install GoogleTest in a way that `pkg-config` knows
where it is. For example, on MacOS just `brew install googletest`.
Then run `make test`.
