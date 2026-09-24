# Describe a dataset's spatial image + coordinates for `scroll_build()`

Names the tissue image (Seurat FOV / `@images` slot) whose spot
coordinates and H&E array should be exported. Pass the result as
`scroll_build(..., spatial = spatial_spec(...))` to add a `spatial`
embedding (from `GetTissueCoordinates()`), bake the tissue image, and
enable the Spatial panel. Defaults auto-detect the first image, so a
standard Visium object needs only `spatial_spec()`.

## Usage

``` r
spatial_spec(image = NULL, name = "spatial", image_asset = TRUE)
```

## Arguments

- image:

  Name of the image / FOV to export (an entry of
  `SeuratObject::Images(object)`). `NULL` (default) uses the first one.

- name:

  Embedding name to give the exported coordinates (default `"spatial"`).
  Appears as a reduction in DimPlot/FeaturePlot and is the Spatial
  panel's coordinate system.

- image_asset:

  Whether to bake the tissue image as an underlay asset (default
  `TRUE`). `FALSE` exports coordinates only (the Spatial panel then
  draws points on a blank background).

## Value

A `scroll_spatial_spec` list.

## See also

[`scroll_build()`](https://shihanli92.github.io/scroll/reference/scroll_build.md)
