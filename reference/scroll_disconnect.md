# Close a scroll project handle

Drops the handle's cached datasets. There is no live database connection
to close (queries read Parquet directly via arrow), so this is a safe
no-op that may be called more than once.

## Usage

``` r
scroll_disconnect(con)
```

## Arguments

- con:

  A handle from
  [`scroll_connect()`](https://shihanli92.github.io/scroll/reference/scroll_connect.md).

## Value

`invisible(NULL)`.
