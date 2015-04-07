

# Introduction #

Gerardus provides some Matlab functions to load, process and save images in [SCI NRRD](http://teem.sourceforge.net/nrrd/format.html) format. This is the native format used by the platform [Seg3D](http://www.sci.utah.edu/cibc/software/42-seg3d.html).

This page contains some notes on the SCI NRRD format.


# Details #

## Rows and columns correspondence to x, y-coordinates ##

**Summary:** The only thing that changes is that we follow the Matlab convention of having the y-coordinate as rows, and the x-coordinate as columns. The NRRD parameters are ordered so that they are consistent with the image data.

### Seg3D's convention for voxel coordinates ###

The way Seg3D stores images in .mat files establishes the following correspondence:

> x <=> rows

> y <=> columns

> z <=> slice

### Matlab and Gerardus' convention for voxel coordinates ###

Note that the Seg3D convention is different from the Matlab convention of using columns for the x-coordinate and rows for the y-coordinate.

Because Matlab is so widely used, in order to make this discrepancy invisible to the user, functions `scinrrd_load()` and `scinrrd_save()` swap the rows and columns accordingly.

That is, data loaded using `scinrrd_load()` follows the convention

> y <=> rows

> x <=> columns

> z <=> slice

The `nrrd.axis` parameters follow the same convention. For example,

> nrrd.axis.spacing(1) <=> dy

> nrrd.axis.spacing(2) <=> dx

> nrrd.axis.spacing(3) <=> dz

This way, the NRRD parameters are consistent with the image data.

This convention is assumed by all the other [Gerardus Matlab](http://code.google.com/p/gerardus/source/browse/#svn/trunk/matlab) functions.

### Everyone's convention for point coordinates ###

Point coordinates are always given as

> (x, y, z)

in Matlab, Seg3D, etc.

For example, if you want to plot a curve on an image in Matlab, you use

```
imagesc(im)
hold on
plot(x, y)
```

where the pixels in `im` have implicitly the coordinates as (y, x), while the curve is plotted as (x, y).

That's why Gerardus functions `scinrrd_index2world()` and `scinrrd_world2index()` by default consider _index_ as (y, x), and _world_ as (x, y).

## Meaning of Offset, point coordinate and index values ##

Correspondence between header fields in the MetaImage (`.mha`) file and in the Matlab `nrrd` struct.

| **MetaImage (`.mha`)** | **Matlab `nrrd` struct** | Comment |
|:-----------------------|:-------------------------|:--------|
| Offset                | `nrrd.axis.min`        | Left edge of first voxel |
| --                    | `nrrd.axis.max`        | Left edge of last voxel |
| ElementSpacing       | `nrrd.axis.spacing`    | Voxel size |
| DimSize               | `nrrd.axis.size`       | Number of voxels |

Note that `nrrd.axis.min` and `nrrd.axis.max` both indicate the "left" side of the first and last voxels.

For example, suppose we have 3 voxels:

| V1 | V2  | V3 |
|:---|:----|:---|

Let `nrrd.axis.min` = `[0 2 1]`, and `nrrd.axis.spacing` = `[3 3 3]`.

Then, the first voxel is defined by the corners
  * `[0 2 1]` = `nrrd.axis.min`
  * `[3 5 4]` = `nrrd.axis.min` + `nrrd.axis.spacing`

The last voxel is defined by the corners
  * `[6 8 7]` = `nrrd.axis.max` = `nrrd.axis.min` + `(nrrd.axis.size - 1) * nrrd.axis.spacing`
  * `[9 11 10]` = `nrrd.axis.min` + `nrrd.axis.size * nrrd.axis.spacing`


Note that the centre of the voxel is 0.5 voxels away from the vertex.

**Point coordinates** are continuous over the whole image volume, while **Indices** are discrete in Seg3D.

In Gerardus, functions `scinrrd_index2world()` and `scinrrd_world2index()` don't round index values, to avoid introducing numerical errors. To obtain true index values, use Matlab's function `round()`.