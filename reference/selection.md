# Variable selection from a stability path

S3 generic for variable selection based on a
[StabilityPath](https://jchiquet.github.io/quadrupen/reference/StabilityPath.md)
object, as introduced by Meinshausen and Buhlmann (2010).

## Usage

``` r
selection(object, ...)

# S3 method for class 'StabilityPath'
selection(
  object,
  sel_mode = c("rank", "PFER"),
  cutoff = 0.75,
  PFER = 2,
  nvarsel = NULL,
  ...
)
```

## Arguments

- object:

  a
  [StabilityPath](https://jchiquet.github.io/quadrupen/reference/StabilityPath.md)
  object.

- ...:

  not used, only here for S3 compatibility.

- sel_mode:

  a character string, either `"rank"` or `"PFER"`. Default is `"rank"`.

- cutoff:

  probability threshold for `sel_mode = "PFER"`. Default is `0.75`.

- PFER:

  per-family error rate to control for `sel_mode = "PFER"`. Default is
  `2`.

- nvarsel:

  number of variables to select for `sel_mode = "rank"`. Default is
  `floor(nobs / log(nvar))`.

## Value

an integer vector of selected variable indices.

## Methods (by class)

- `selection(StabilityPath)`: S3 method for variable selection from a
  [StabilityPath](https://jchiquet.github.io/quadrupen/reference/StabilityPath.md)

## References

N. Meinshausen and P. Buhlmann (2010). Stability Selection, JRSS(B).

## See also

[`stability()`](https://jchiquet.github.io/quadrupen/reference/stability.md)
