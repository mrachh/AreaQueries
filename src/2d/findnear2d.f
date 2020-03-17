c              this file contains the near field routines 
c              to be used for quadrature generation
c
c              findnear2d - find near field of a collection
c                         of sources with extent and targets
c 
c              findnearmem2d - memory management routine for find
c                            near
c
c              O(n^2) analogs of findnear and findnearmem
c
c              findnearslow2d - find near field of a collection
c                         of sources with extent and targets
c 
c              findnearslowmem2d - memory management routine for find
c                            near
c
c
c-------------------------------------------------------------------
      subroutine findnearslow2d(xyzs,ns,rads,targets,nt,row_ptr,
     1   col_ind)
c
cc      identify targets which are contained in 
c       |xyzs(:,i)-targets(:,j)|<=rads(i)
c       This is a dense O(ns \times nt) algorithm
c
c       Calling sequence variables
c       xyzs    in: real *8 (2,ns)
c               location of the sources
c
c       ns      in: integer
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       targets in: real *8(2,nt)
c               location of the targets
c
c       nt      in: integer
c               number of targets
c
c
c       OUTPUT
c       row_ptr    out: integer(nt+1)
c                rowptr(i) is the starting point in the iflg
c                array for the list of sources
c                relevant for target i
c
c       col_ind  out: integer(nnz) (nnz is computed using mem routine)
c                col_ind(row_ptr(i):row_ptr(i+1)) is the list
c                of sources which satisfy
c                d(s_{j},t_{i}) <= rs_{j}
c               
c-------------------------------

      implicit real *8 (a-h,o-z)
      real *8 xyzs(2,ns),targets(2,nt),rads(ns)
      integer row_ptr(*),col_ind(*)
      integer, allocatable :: nlst(:)

      allocate(nlst(nt))

      do i=1,nt
        nlst(i) = 0
        do j=1,ns
          rr = (xyzs(1,j)-targets(1,i))**2 + 
     1          (xyzs(2,j)-targets(2,i))**2 
          if(rr.le.rads(j)**2) then
            nlst(i) = nlst(i) + 1
          endif
        enddo
      enddo

      row_ptr(1) = 1
      do i=1,nt
        row_ptr(i+1) = row_ptr(i)+nlst(i)
      enddo

      do i=1,nt
        nlst(i) = 0
      enddo

      do i=1,nt
        do j=1,ns
          rr = (xyzs(1,j)-targets(1,i))**2 + 
     1         (xyzs(2,j)-targets(2,i))**2 

      
          if(rr.le.rads(j)**2) then
            col_ind(row_ptr(i)+nlst(i)) = j
            nlst(i) = nlst(i) + 1
          endif
        enddo
      enddo



      return
      end
c
c
c
c
c
c------------------------------------------------------------------      

      subroutine findnearslowmem2d(xyzs,ns,rads,targets,nt,nnz)
c
cc      identify all sources which are contained in 
c       |xyzs(:,i)-targets(:,j)|<=rads(i)
c       This is a dense O(ns \times nt) algorithm.
c
c       Calling sequence variables
c       xyzs    in: real *8 (2,ns)
c               location of the sources
c
c       ns      in: integer
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       targets in: real *8(2,nt)
c               location of the targets
c
c       nt      in: integer
c               number of targets
c
c
c       OUTPUT
c       nnz     out: integer
c               number of elements in the flag array
c-------------------------------

      implicit real *8 (a-h,o-z)
      real *8 xyzs(2,*),targets(2,*),rads(*)
      integer nnz

      nnz = 0

      do i=1,nt
        do j=1,ns
          rr = (xyzs(1,j)-targets(1,i))**2 + 
     1         (xyzs(2,j)-targets(2,i))**2 

          if(rr.le.rads(j)**2) then
            nnz = nnz+1
          endif
        enddo
      enddo

      return
      end
c
c
c
c
c
c------------------------------------------------------------------      
      subroutine findnearmem2d(xyzs,ns,rads,targets,nt,nnz)
c
cc      identify all sources which are contained in 
c       |xyzs(:,i)-targets(:,j)|<=rads(i).
c       We use hung lists to identify a nearly minimal set of
c       sources to loop over relevant for each target.
c
c       This is a memory management routine for findnearmem
c       
c
c       Calling sequence variables
c
c       xyzs    in: real *8 (2,ns)
c               location of the sources
c
c       ns      in: integer
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       targets in: real *8(2,nt)
c               location of the targets
c
c       nt      in: integer
c               number of targets
c
c       radt    in: real *8(ns)
c               radii associated with the targets
c
c       targets in: real *8(2,nt)
c               location of the targets
c
c       nt      in: integer
c               number of targets
c
c       OUTPUT
c       nnz     out: integer
c               number of elements in the flag array
c
c               
c-------------------------------

       implicit real *8 (a-h,o-z)
       real *8 xyzs(2,ns),rads(ns),targets(2,nt)
       real *8 xtmp(2)
       real *8, allocatable :: rstmp(:)

       integer nnz

cc       tree variables
c 
       real *8, allocatable :: centers(:,:),boxsize(:)
       integer, allocatable :: ilevel(:)
       integer, allocatable :: itree(:)
       integer ipointer(30),ltree

       logical res



       allocate(rstmp(ns))



       do i=1,ns
         rstmp(i) = 2*rads(i)
       enddo


       rttmp = 0


       idivflag = 0
       ndiv = 2
       isep = 1
       nlmax = 200
       nbmax = 0

       nlevels = 0
       nboxes = 0
       mhung = 0
       mnbors = 0
       mnlist1 = 0
       mnlist2 = 0
       mnlist3 = 0
       mnlist4 = 0
       ntmp = 0


       call maketree2dmem(ier,xyzs,ns,rstmp,targets,nt,xtmp,ntmp,rttmp,
     1    idivflag,ndiv,nlmax,nbmax,nlevels,nboxes,mhung,ltree)


       allocate(centers(2,nboxes),itree(ltree),boxsize(0:nlevels))
  
       call maketree2d(xyzs,ns,rstmp,targets,nt,xtmp,ntmp,rttmp,
     1    idivflag,ndiv,mhung,nlevels,nboxes,centers,boxsize,itree,
     2    ltree,ipointer,mnlist1,mnlist2,mnlist3,mnlist4)


       allocate(ilevel(nboxes))

       do ilev=0,nlevels
         do ibox=itree(2*ilev+1),itree(2*ilev+2)
           ilevel(ibox) = ilev
         enddo
       enddo

       nnz = 0


       do ilev=0,nlevels
         do ibox=itree(2*ilev+1),itree(2*ilev+2)
           nchild = itree(ipointer(3)+ibox-1)
           if(nchild.eq.0) then 
             itstart = itree(ipointer(12)+ibox-1)
             itend = itree(ipointer(13)+ibox-1)
             do itt = itstart,itend
               itarg = itree(ipointer(6)+itt-1)
               nhunglistsrc = itree(ipointer(28)+ibox-1)
               do ii=1,nhunglistsrc
                 iss = itree(ipointer(29)+mhung*(ibox-1)+ii-1)
                 rr = (targets(1,itarg)-xyzs(1,iss))**2+ 
     1                 (targets(2,itarg)-xyzs(2,iss))**2
                 if(rr.le.rads(iss)**2) nnz = nnz + 1              
               enddo

               nlist1 = itree(ipointer(18)+ibox-1)
               do j=1,nlist1
                 jbox = itree(ipointer(19)+mnlist1*(ibox-1)+j-1)
                 isstart = itree(ipointer(10)+jbox-1)
                 isend = itree(ipointer(11)+jbox-1)
                 do ii = isstart,isend
                   iss = itree(ipointer(5)+ii-1)
                   rr = (targets(1,itarg)-xyzs(1,iss))**2+ 
     1                    (targets(2,itarg)-xyzs(2,iss))**2
                   if(rr.le.rads(iss)**2) nnz = nnz + 1
                 enddo
               enddo
             enddo
           endif
         enddo
       enddo


       return
       end
c
c
c
c
c
c-----------------------------------------------------
      subroutine findnear2d(xyzs,ns,rads,targets,nt,row_ptr,
     1       col_ind) 
c     
cc      identify all sources which are contained in 
c       |xyzs(:,i)-targets(:,j)|<=rads(i).
c       We use hung lists to identify a nearly minimal set of
c       sources to loop over relevant for each target.
c       
c
c       Calling sequence variables
c
c       xyzs    in: real *8 (2,ns)
c               location of the sources
c
c       ns      in: integer
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       targets in: real *8(2,nt)
c               location of the targets
c
c       nt      in: integer
c               number of targets
c
c       row_ptr    out: integer(nt+1)
c                row_ptr(i) is the starting point in the col_ind
c                array for the list of sources
c                relevant for target i
c
c       col_ind     out: integer(nnz) (nnz is computed using
c                                      findnearmem)
c                col_ind(row_ptr(i):row_ptr(i+1)-1) is the list
c                of sources which satisfy
c                d(s_{j},t_{i}) <= rs_{j}
c               
c               
c-------------------------------
       implicit real *8 (a-h,o-z)
       real *8 xyzs(2,*),rads(*),targets(2,*)
       real *8, allocatable :: rstmp(:)

       integer row_ptr(*),col_ind(*)
       integer, allocatable :: nlst(:)
c
cc       tree variables
c 
       real *8, allocatable :: centers(:,:),boxsize(:)
       integer, allocatable :: itree(:)
       integer ipointer(30), ltree
       integer, allocatable :: ilevel(:)

       allocate(rstmp(ns))

       do i=1,ns
         rstmp(i) = 2*rads(i)
       enddo

       rttmp = 0.0d0



       idivflag = 0
       ndiv = 2
       isep = 1
       nlmax = 200
       nbmax = 0

       nlevels = 0
       nboxes = 0
       mhung = 0
       ntmp = 0
       mnlist1 = 0
       mnlist2 = 0
       mnlist3 = 0
       mnlist4 = 0
       mnbors = 0


       call maketree2dmem(ier,xyzs,ns,rstmp,targets,nt,xtmp,ntmp,rttmp,
     1    idivflag,ndiv,nlmax,nbmax,nlevels,nboxes,mhung,ltree)


       allocate(centers(2,nboxes),itree(ltree),boxsize(0:nlevels))
  
       call maketree2d(xyzs,ns,rstmp,targets,nt,xtmp,ntmp,rttmp,
     1    idivflag,ndiv,mhung,nlevels,nboxes,centers,boxsize,itree,
     2    ltree,ipointer,mnlist1,mnlist2,mnlist3,mnlist4)

       allocate(ilevel(nboxes))

       do ilev=0,nlevels
          do ibox=itree(2*ilev+1),itree(2*ilev+2)
              ilevel(ibox) = ilev
          enddo
       enddo

       allocate(nlst(nt))


       do i=1,nt
          nlst(i) = 0
       enddo

       do ilev=0,nlevels
         do ibox=itree(2*ilev+1),itree(2*ilev+2)
           nchild = itree(ipointer(3)+ibox-1)
           if(nchild.eq.0) then
             itstart = itree(ipointer(12)+ibox-1)
             itend = itree(ipointer(13)+ibox-1)

             do it = itstart,itend
               itarg = itree(ipointer(6)+it-1)

               nhunglistsrc = itree(ipointer(28)+ibox-1)
               do ii=1,nhunglistsrc
                 is = itree(ipointer(29)+mhung*(ibox-1)+ii-1)
                 rr = (targets(1,itarg)-xyzs(1,is))**2+ 
     1                 (targets(2,itarg)-xyzs(2,is))**2
                 if(rr.le.rads(is)**2) 
     1               nlst(itarg) = nlst(itarg) + 1              
               enddo

               nlist1 = itree(ipointer(18)+ibox-1)
               do j=1,nlist1
                 jbox = itree(ipointer(19)+mnlist1*(ibox-1)+j-1)
                 isstart = itree(ipointer(10)+jbox-1)
                 isend = itree(ipointer(11)+jbox-1)
                 do ii = isstart,isend
                   is = itree(ipointer(5)+ii-1)
                   rr = (targets(1,itarg)-xyzs(1,is))**2+ 
     1                     (targets(2,itarg)-xyzs(2,is))**2
                   if(rr.le.rads(is)**2) 
     1                    nlst(itarg) = nlst(itarg) + 1
                 enddo
               enddo
             enddo
           endif
         enddo
       enddo

       row_ptr(1) = 1
       do i=1,nt
          row_ptr(i+1) = row_ptr(i) + nlst(i)
          nlst(i) = 0
       enddo



       do ilev=0,nlevels
         do ibox=itree(2*ilev+1),itree(2*ilev+2)
           nchild = itree(ipointer(3)+ibox-1)
           if(nchild.eq.0) then
             itstart = itree(ipointer(12)+ibox-1)
             itend = itree(ipointer(13)+ibox-1)

             do it = itstart,itend
               itarg = itree(ipointer(6)+it-1)
               nhunglistsrc = itree(ipointer(28)+ibox-1)
               do ii=1,nhunglistsrc
                 is = itree(ipointer(29)+mhung*(ibox-1)+ii-1)
                 rr = (targets(1,itarg)-xyzs(1,is))**2+ 
     1                 (targets(2,itarg)-xyzs(2,is))**2
                 if(rr.le.rads(is)**2) then
                   col_ind(row_ptr(itarg)+nlst(itarg)) =is
                   nlst(itarg) = nlst(itarg)+1
                 endif
               enddo

               nlist1 = itree(ipointer(18)+ibox-1)
               do j=1,nlist1
                 jbox = itree(ipointer(19)+mnlist1*(ibox-1)+j-1)
                 isstart = itree(ipointer(10)+jbox-1)
                 isend = itree(ipointer(11)+jbox-1)
                 do ii = isstart,isend
                   is = itree(ipointer(5)+ii-1)
                   rr = (targets(1,itarg)-xyzs(1,is))**2+ 
     1                  (targets(2,itarg)-xyzs(2,is))**2
                   if(rr.le.rads(is)**2) then 
                     col_ind(row_ptr(itarg)+nlst(itarg)) =is
                     nlst(itarg) = nlst(itarg)+1
                   endif
                 enddo
               enddo
             enddo
           endif
         enddo
       enddo


       return
       end
c
c
