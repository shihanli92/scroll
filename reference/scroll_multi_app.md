# Build a multi-dataset scroll explorer

Mounts several built scroll projects behind one app: a tab per dataset,
each with its own app bar, manifest-driven panels, and query handle.
Every dataset's UI is namespaced (via
[`shiny::NS()`](https://rdrr.io/pkg/shiny/man/NS.html)) so the projects
coexist without id collisions; the built-in and custom panels are the
same registry for all.

## Usage

``` r
scroll_multi_app(projects)
```

## Arguments

- projects:

  A named list/vector of built project directories. Names are used as
  tab labels (falling back to each project's `config.yaml` title, then a
  generic label).

## Value

A `shiny.appobj`.

## Details

Custom panels that read baked assets should resolve them from the
per-dataset handle (`data$dir`) rather than a global option, so each tab
reads its own project's files.

The `config.yaml` panel override (`panels:` / `exclude_panels:`) is
honoured by
[`scroll_app()`](https://shihanli92.github.io/scroll/reference/scroll_app.md);
here the tab panel list is shared, so the panels are the gated union
across all mounted projects.
