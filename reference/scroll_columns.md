# Metadata columns of a given type

Convenience accessor over a scroll data handle's manifest, for
populating custom-panel controls.

## Usage

``` r
scroll_columns(data, type = c("categorical", "numeric", "any"), view = NULL)
```

## Arguments

- data:

  A scroll data handle (as passed to a panel `ui`/`server`).

- type:

  `"categorical"`, `"numeric"`, or `"any"`.

- view:

  Optional subset-view name; when given, subset-scoped columns are
  restricted to that view (global columns always appear). `NULL`
  (default) returns the whole-dataset column set.

## Value

A character vector of column names.
