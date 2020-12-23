# AreaQueries

Author: Manas Rachh, Flatiron Institute.

This repository contains Fortran90 routines for fast identification of
the list of test (target) points contained within each of a given
collection of (source) balls, in 2D or 3D. The results are returned in
compressed sparse row (CSR) format, ie, indices of balls that contain
each target. The algorithm used a hung pruned level-restricted
adaptive quad-tree in 2D (or oct-tree in 3D), with OpenMP parallelization.

### Example

We present the 2D case (3D is analogous). Given inputs:

  `s` : 2-by-`ns` array giving centers of source balls in R<sup>2</sup>  
  `r` : length-`ns` list of radii of source balls  
  `t` : 2-by-`nt` array giving coords of test points in R<sup>2</sup>  

Then calling:
```
      call findnearmem2d(s,ns,r,t,nt,nnz)
      allocate(col_ind(nnz))
      call findnear2d(s,ns,r,t,nt,row_ptr,col_ind)
```
Writes out the following:

  `nnz` (written by `findnearmem2d`) : total number of target-ball pairs found  
  `row_ptr` : length-`(nt+1)` list of start indices in `col_ind` array for each target  
  `col_ind` : length-`nnz` list of source ball indices intersecting each target  

This may be interpreted as the CSR storage of a `nt`-by-`ns`
sparse matrix with `nnz` nonzero entries.
The element `(i,j)` is nonzero iff target `i` is within ball `j`, ie,
within distance `r(j)` of source `j`.

The task performed may be also understood by reading the naive codes
`findnearslowmem2d` and `findnearslow2d`.


### Contents of package

`src/2d` : 2D source code  
`src/3d` : 3D source code  
`src/common` : sparse matrix format conversions, low-level utilities including by other authors  
`test` : test drivers and makefiles  

### Testing

```
cd test
./testall.sh
```
should compile and run `test2d` and `test3d`. Each reports timing tests vs the naive method (taking around 1 sec to run), and successful completion of test of correctness of indices found.
