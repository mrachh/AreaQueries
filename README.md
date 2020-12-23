# AreaQueries

Author: Manas Rachh, Flatiron Institute.

This repository contains Fortran90 routines for fast identification of list of target points contained within specified neighborhoods of a given collection of source points. The results are returned in compressed sparse row (CSR) format. There are 2D and 3D versions.

### 2D example

Given:
  `s` is 2-by-`ns` array giving coords of points in R2
  ``r`` is length ``ns`` list of radii corresponding to above points  
  ``t`` is 2-by-``nt`` array giving coords of points in R2  

Then:
```
      call findnearmem2d(s,ns,r,t,nt,nnz)
      allocate(col_ind(nnz))
      call findnear2d(src,ns,rads,targ,nt,row_ptr,col_ind)
```



### Contents of package

``src/2d`` : 2D source code
``src/3d`` : 3D source code
``src/common`` : sparse matrix format conversions, utilities by other authors
``test`` : test drivers and makefiles

### Testing

