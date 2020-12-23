# AreaQueries

Author: Manas Rachh, Flatiron Institute.

This repository contains Fortran90 routines for fast identification of list of target points contained within specified neighborhoods of a given collection of source points. The results are returned in compressed sparse row (CSR) format. There are 2D and 3D versions.

### 2D example

Given inputs:

  `s` : 2-by-`ns` array giving coords of source points in R2  
  `r` : length `ns` list of source radii  
  `t` : 2-by-`nt` array giving coords of target points in R2  

Then calling:
```
      call findnearmem2d(s,ns,r,t,nt,nnz)
      allocate(col_ind(nnz))
      call findnear2d(src,ns,rads,targ,nt,row_ptr,col_ind)
```
Writes out the following:

  `nnz` (written by `findnearmem2d`) : total number of within-radius source-target pairs found  
  `row_ptr` : length `nt+1` list of start indices in `col_ind` array for each target  
  `col_ind` : length `nnz` list of source indices within their respective radii of each target  

This may be interpreted as the CSR storage of a `nt`-by-`ns`
sparse matrix with `nnz` nonzero entries.
The element `(i,j)` is nonzero when target 'i' is within radius
`r(j)` of source `j`.

### Contents of package

``src/2d`` : 2D source code
``src/3d`` : 3D source code
``src/common`` : sparse matrix format conversions, utilities by other authors
``test`` : test drivers and makefiles

### Testing

