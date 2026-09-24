# Open a scroll project for querying

Returns a lightweight handle over a built project's Parquet store. The
handle caches one
[`arrow::Dataset`](https://arrow.apache.org/docs/r/reference/Dataset.html)
per assay (built lazily on first use), so feature lookups avoid
re-scanning the partition directory. Release it with
[`scroll_disconnect()`](https://shihanli92.github.io/scroll/reference/scroll_disconnect.md).

## Usage

``` r
scroll_connect(dir)
```

## Arguments

- dir:

  A scroll project directory (must contain `expr/`).

## Value

A `scroll_con` handle. Close it with
[`scroll_disconnect()`](https://shihanli92.github.io/scroll/reference/scroll_disconnect.md).
