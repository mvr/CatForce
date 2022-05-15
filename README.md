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

| Line                  | Parameter             | Description                                                            |
|-----------------------|-----------------------|------------------------------------------------------------------------|
| `max-gen`             | `n`                   | Overall maximum generation                                             |
| `start-gen`           | `n`                   | The min gen the first encounter must happen by                         |
| `last-gen`            | `n`                   | The max gen the _first_ encounter is allowed                           |
| `num-catalyst`        | `n`                   | The number of catalyst to place                                        |
| `stable-interval`     | `n`                   | Gens the catalysts must remain untouched to be considered stable       |
| `search-area`         | `x y w h`             | Search area                                                            |
| `pat`                 | `rle`                 | The active pattern                                                     |
|                       | `(dx dy)`             | Optional offset                                                        |
| `cat `                | `rle`                 | A catalyst                                                             |
|                       | `max-active`          | Number of generations in a row the catalyst may be missing             |
|                       | `dx dy`               | Offset to centre the catalyst (typically negative)                     |
|                       | `symmetries-char`     | Character denoting the symmetry of the catalyst                        |
|                       | `(forbidden rle x y)` | Optional forbidden pattern around the catalyst                         |
| `output`              | `filename`            | Output filename                                                        |
| `filter`              | `gen rle dx dy`       | Filter that must be matched for the solution to be accepted            |
| `filter`              | `min-max rle dx dy`   | Filter in a range of generations.                                      |
| `full-report`         | `filename`            | Output filename for solutions ignoring all pattern filters             |
| `max-category-size`   | `n`                   | Maximum output row length before more solutions are dropped            |
| `fit-in-width-height` | `w h`                 | Only allow solutions where all catalysts fit in a `w` by `h` rectangle |
| `combine-results`     |                       | See below                                                              |
| `symmetry`            |                       | See below                                                              |



**Encounter**: An active cell from the input pattern is present in
the immediate neighbourhood of the catalyst.

**Catalyst Symmetry**: A character specifying what symmetries are
applied to the catalyst:
- `|` mirror by y.
- `-` mirror by x.
- `+` mirror by x, y and both. 
- `/` diagonal mirror.
- `x` make for 180 degree symmetrical catalysts.
- `*` all 8 transformations.

**Forbidden**: Will search for `rle` in the same location as the
catalyst. If `rle` is matched in any generation, the solution is
omitted from the report. You can see forbidden as "catalyst attached
filter", i.e. the forbidden pattern is based on catalyst location.
Made to exclude eaters/boat-bits and similar unwanted "garbage". You
may have several forbidden patterns per catalyst. See `3.in` for
an example.

**Filters**: one can use several filters. Every filter will be
checked if successful catalyst was found. Each filter will assume live
cells in `rle` and dead cell in close proximity neighbourhood to the
pattern in `rle`.

Combining Results
---

`combine-results yes [<survive-0> <survive-1> ...]`

If this feature is enabled the search will at first ignore all filters
and survival inputs, and will search all the possible catalysts. Then
it will try to combine all the found catalysts in all possible
combinations, and only then will filter by `survive-i` and apply the
filters to exclude them from the final report.
 
This feature will generate report as follows:
 
- `output.rle` - all the possible catalysts.
- `output.rle_Combined*.rle` - will generate all combined reports.
- `output.rle_Final.rle` - the final report. **This is the main output.**
 
Optional survival filter per "iteration" are added. Combine works as
follows: each time it start from the initial search results (combine
by default uses survive count = 1), and tries to add catalyst from
those results. Sometimes one could get explosion, if the interaction
is very potent. So filter is added to limit the combine, by surviving
count (if something doesn't survive with two catalyst for 5
iterations, it's probably junk - so CatForce will filter it on the
second combine iteration and not in the end).
 
This allows faster and more efficient combine operation with very
potent conduits which otherwise would overflow the system, with many
useless catalysts.
 
**NOTE** Recommended for use only for `num-catalyst` = 1 or 2

**NOTE** See 4.in file for example. 
 
**NOTE** CatForce will use the last `survive-i` as the default from
that point on. If you don't enter any numbers it will use survival
count 1, and will filter only when finish all possible combinations.

Symmetric Searches
---

`symmetry s`

Automatically apply the symmetry `s` to the active pattern and
catalysts. The lines of symmetry are always `x=0` and `y=0`, and
so if you want a different offset you will have to change the `dx dy`
parameters of the active pattern.

Options are:
- `horizontal`
- `horizontaleven`
- `diagonal`
- `diagonalevenx` 
- `diagonalevenboth`
 
