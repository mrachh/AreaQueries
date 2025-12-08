c              this file contains the near field routines 
c              to be used for quadrature generation
c
c              findnear3d - find near field of a collection
c                         of sources with extent and targets
c 
c              findnearmem3d - memory management routine for find
c                            near
c
c              O(n^2) analogs of findnear and findnearmem
c
c              findnearslow3d - find near field of a collection
c                         of sources with extent and targets
c 
c              findnearslowmem3d - memory management routine for find
c                            near
c
c
c-------------------------------------------------------------------
      subroutine findnearslow3d(xyzs,ns,rads,ndt,targets,nt,row_ptr,
     1  col_ind)
c
cc      identify targets which are contained in 
c       |xyzs(:,i)-targets(:,j)|<=rads(i)
c       This is a dense O(ns \times nt) algorithm
c
c       Calling sequence variables
c       xyzs    in: real *8 (3,ns)
c               location of the sources
c
c       ns      in: integer
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       ndt     in: integer
c               leading dimension of targets array
c
c       targets in: real *8(ndt,nt)
c               target info.
c               The first three components must be xyz coordinates
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
      real *8 xyzs(3,ns),targets(ndt,nt),rads(ns)
      integer row_ptr(*),col_ind(*)
      integer, allocatable :: nlst(:)

      allocate(nlst(nt))

      do i=1,nt
        nlst(i) = 0
        do j=1,ns
          rr = (xyzs(1,j)-targets(1,i))**2 + 
     1          (xyzs(2,j)-targets(2,i))**2 +
     2          (xyzs(3,j)-targets(3,i))**2
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
     1         (xyzs(2,j)-targets(2,i))**2 +
     2        (xyzs(3,j)-targets(3,i))**2

      
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

      subroutine findnearslowmem3d(xyzs,ns,rads,ndt,targets,nt,nnz)
c
cc      identify all sources which are contained in 
c       |xyzs(:,i)-targets(:,j)|<=rads(i)
c       This is a dense O(ns \times nt) algorithm.
c
c       Calling sequence variables
c       xyzs    in: real *8 (3,ns)
c               location of the sources
c
c       ns      in: integer
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c      
c       ndt     in: integer
c               leading dimension of target array
c
c       targets in: real *8(ndt,nt)
c               target info
c               the first three components must be
c               the xyz coordinates
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
      real *8 xyzs(3,*),targets(ndt,*),rads(*)
      integer nnz

      nnz = 0

      do i=1,nt
        do j=1,ns
          rr = (xyzs(1,j)-targets(1,i))**2 + 
     1         (xyzs(2,j)-targets(2,i))**2 +
     2         (xyzs(3,j)-targets(3,i))**2

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
      subroutine findnearmem3d(xyzs,ns,rads,ndt,targets,nt,nnz)
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
c       xyzs    in: real *8 (3,ns)
c               location of the sources
c
c       ns      in: integer
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       ndt     in: integer
c               leading dimension of target array
c
c       targets in: real *8(ndt,nt)
c               target info
c               the first three components must be
c               the xyz coordinates
c
c       nt      in: integer
c               number of targets
c
c       radt    in: real *8(ns)
c               radii associated with the targets
c
c       targets in: real *8(3,nt)
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
       real *8 xyzs(3,*),rads(*),targets(ndt,*)

       idivflag = 0
       ndiv = 2
       call findnearmem3d_tparams(xyzs,ns,rads,ndt,targets,nt,idivflag,
     1    ndiv,nnz)

       return
       end
c
c
c
c
c
c-----------------------------------------------------
      subroutine findnear3d(xyzs,ns,rads,ndt,targets,nt,row_ptr,
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
c       xyzs    in: real *8 (3,ns)
c               location of the sources
c
c       ns      in: integer
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       ndt     in: integer
c               leading dimension of target array
c
c       targets in: real *8(ndt,nt)
c               target info
c               the first three components must be
c               the xyz coordinates
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
       real *8 xyzs(3,*),rads(*),targets(ndt,*)
       integer row_ptr(*),col_ind(*)

       idivflag = 0
       ndiv = 2

       call findnear3d_tparams(xyzs, ns, rads, ndt, targets, 
     1    nt, idivflag, ndiv, row_ptr, col_ind)


      return
      end
c
c------------------------------------------------------------------      
      subroutine findnearmem3d_tparams(xyzs,ns,rads,ndt,targets,
     1    nt, idivflag, ndiv, nnz)
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
c       xyzs    in: real *8 (3,ns)
c               location of the sources
c
c       ns      in: integer
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       ndt     in: integer
c               leading dimension of target array
c
c       targets in: real *8(ndt,nt)
c               target info
c               the first three components must be
c               the xyz coordinates
c
c       nt      in: integer
c               number of targets
c
c       radt    in: real *8(ns)
c               radii associated with the targets
c
c       targets in: real *8(3,nt)
c               location of the targets
c
c       nt      in: integer
c               number of targets
c
c       idivflag in:integer
c        subdivision criterion
c          * divflag = 0 -> subdivide on sources only
c          * idivflag = 1 -> subdivide on targets only
c          * idivflag = 2 -> subdivide on max(sources+targets)
c       ndiv    in: integer
c                 max number of points in leaf box
c       OUTPUT
c       nnz     out: integer
c               number of elements in the flag array
c
c               
c-------------------------------

       implicit real *8 (a-h,o-z)
       real *8 xyzs(3,*),rads(*),targets(ndt,*)
       real *8, allocatable :: rstmp(:)

       integer, allocatable :: nlst(:)
c
cc       tree variables
c 
       real *8, allocatable :: centers(:,:),boxsize(:),targs(:,:)
       integer, allocatable :: itree(:)
       integer *8 ipointer(8), ltree
       integer, allocatable :: ilevel(:)
       integer, allocatable :: nlist1(:),list1(:,:)
       integer, allocatable :: nlist2(:),list2(:,:)
       integer, allocatable :: nlist3(:),list3(:,:)
       integer, allocatable :: nlist4(:),list4(:,:)
       integer mnlist1,mnlist2,mnlist3,mnlist4
       integer, allocatable :: isrcper(:),isrcse(:,:)
       integer, allocatable :: itargper(:),itargse(:,:)

       allocate(rstmp(ns))

C$OMP PARALLEL DO DEFAULT(SHARED)
       do i=1,ns
         rstmp(i) = 2*rads(i)
       enddo
C$OMP END PARALLEL DO       
       allocate(targs(3,nt))

C$OMP PARALLEL DO DEFAULT(SHARED)       
       do i=1,nt
         targs(1,i) = targets(1,i)
         targs(2,i) = targets(2,i)
         targs(3,i) = targets(3,i)
       enddo
C$OMP END PARALLEL DO       

       isep = 1
       nlmax = 51
       nlmin = 0
       iper = 0
       ifunif = 0
       nbmax = 0

       nlevels = 0
       nboxes = 0
       mnlist1 = 0
       mnlist2 = 0
       mnlist3 = 0
       mnlist4 = 0
       mnbors = 27

       call pts_tree_mem(xyzs,ns,targs,nt,idivflag,ndiv,nlmin,nlmax,
     1  iper,ifunif,nlevels,nboxes,ltree) 

       allocate(centers(3,nboxes),itree(ltree),boxsize(0:nlevels))

       call pts_tree_build(xyzs,ns,targs,nt,idivflag,ndiv,nlmin,
     1   nlmax,iper,ifunif,nlevels,nboxes,ltree,itree,ipointer,centers,
     2   boxsize)
       
       allocate(isrcse(2,nboxes),itargse(2,nboxes))
       allocate(isrcper(ns),itargper(nt))

       call pts_tree_sort(ns,xyzs,itree,ltree,nboxes,nlevels,
     1   ipointer,centers,isrcper,isrcse)
      
       call pts_tree_sort(nt,targs,itree,ltree,nboxes,nlevels,
     1   ipointer,centers,itargper,itargse)

c
c   initialize various tree lists
c
      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0
      mnbors = 27

      isep = 1
      
      call computemnlists(nlevels,nboxes,itree(ipointer(1)),boxsize,
     1  centers,itree(ipointer(3)),itree(ipointer(4)),
     2  itree(ipointer(5)),isep,itree(ipointer(6)),mnbors,
     2  itree(ipointer(7)),iper,mnlist1,mnlist2,mnlist3,mnlist4)
      
      allocate(list1(mnlist1,nboxes),nlist1(nboxes))
      allocate(list2(mnlist2,nboxes),nlist2(nboxes))
      allocate(list3(mnlist3,nboxes),nlist3(nboxes))
      allocate(list4(mnlist4,nboxes),nlist4(nboxes))

      call computelists(nlevels,nboxes,itree(ipointer(1)),boxsize,
     1  centers,itree(ipointer(3)),itree(ipointer(4)),
     2  itree(ipointer(5)),isep,itree(ipointer(6)),mnbors,
     3  itree(ipointer(7)),iper,nlist1,mnlist1,list1,nlist2,
     4  mnlist2,list2,nlist3,mnlist3,list3,
     4  nlist4,mnlist4,list4)

      nnz = 0
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,nchild,ilev,iss)
C$OMP$PRIVATE(isrc,rtest,ilevup,jbox,i,ncoll,l,lbox,xdis,ydis,zdis)
C$OMP$PRIVATE(llev,distest,itt,itarg,rr,jlev) REDUCTION(+:nnz)
C$OMP$SCHEDULE(DYNAMIC)
      do ibox = 1,nboxes
        nchild = itree(ipointer(4)+ibox-1)
        ilev = itree(ipointer(2)+ibox-1)
        if(nchild.eq.0) then
          do iss=isrcse(1,ibox),isrcse(2,ibox)
             isrc = isrcper(iss)
             rtest = rads(isrc)**2
             ilevup = ceiling(log(rstmp(isrc)/boxsize(ilev))/log(2.0d0))
             jbox = ibox
             ilevup = min(ilevup,ilev)

             do i=1,ilevup
               jbox = itree(ipointer(3)+jbox-1)
             enddo
c
c   loop over colleagues
c
             ncoll = itree(ipointer(6)+jbox-1)
             do l=1,ncoll
               lbox = itree(ipointer(7)+mnbors*(jbox-1)+l-1)
c
c  check if lbox does not intersect with region of interest
c 
               xdis = abs(centers(1,lbox)-xyzs(1,isrc))
               ydis = abs(centers(2,lbox)-xyzs(2,isrc))
               zdis = abs(centers(3,lbox)-xyzs(3,isrc))

               llev = itree(ipointer(2)+lbox-1)
               distest = rads(isrc)+boxsize(llev)/2
               if((xdis.le.distest).and.(ydis.le.distest).and.
     1             (zdis.le.distest)) then
                 
                 do itt=itargse(1,lbox),itargse(2,lbox)
                   itarg = itargper(itt)
                   rr = (xyzs(1,isrc)-targs(1,itarg))**2 + 
     1                (xyzs(2,isrc)-targs(2,itarg))**2 + 
     2                (xyzs(3,isrc)-targs(3,itarg))**2
                   if(rr.le.rtest) then
                     nnz = nnz + 1
                   endif
                 enddo
               endif
             enddo

c
c   loop over list1 boxes which are larger
c
             jlev = itree(ipointer(2)+jbox-1)
             do l=1,nlist1(jbox)
               lbox = list1(l,jbox)
               llev = itree(ipointer(2)+lbox-1)
               if(llev.lt.jlev) then
                 do itt=itargse(1,lbox),itargse(2,lbox)
                   itarg = itargper(itt)
                   rr = (xyzs(1,isrc)-targs(1,itarg))**2 + 
     1                (xyzs(2,isrc)-targs(2,itarg))**2 + 
     2                (xyzs(3,isrc)-targs(3,itarg))**2
                   if(rr.le.rtest) then
                     nnz = nnz + 1
                   endif
                 enddo
               endif
             enddo
          enddo
        endif
      enddo
C$OMP END PARALLEL DO      


       return
       end
c
c
c
c
c
c
c
c-----------------------------------------------------
      subroutine findnear3d_tparams(xyzs,ns,rads,ndt,targets,nt,
     1   idivflag,ndiv,row_ptr,col_ind) 
c     
cc      identify all sources which are contained in 
c       |xyzs(:,i)-targets(:,j)|<=rads(i).
c       We use hung lists to identify a nearly minimal set of
c       sources to loop over relevant for each target.
c       
c
c       Calling sequence variables
c
c       xyzs    in: real *8 (3,ns)
c               location of the sources
c
c       ns      in: integer
c               number of sources
c
c       rads    in: real *8(ns)
c               radii associated with the sources
c
c       ndt     in: integer
c               leading dimension of target array
c
c       targets in: real *8(ndt,nt)
c               target info
c               the first three components must be
c               the xyz coordinates
c
c       nt      in: integer
c               number of targets
c       
c       idivflag in:integer
c        subdivision criterion
c          * divflag = 0 -> subdivide on sources only
c          * idivflag = 1 -> subdivide on targets only
c          * idivflag = 2 -> subdivide on max(sources+targets)
c       ndiv    in: integer
c                 max number of points in leaf box
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
       real *8 xyzs(3,*),rads(*),targets(ndt,*)
       real *8, allocatable :: rstmp(:)

       integer row_ptr(*),col_ind(*)
       integer, allocatable :: nlst(:)
c
cc       tree variables
c 
       real *8, allocatable :: centers(:,:),boxsize(:),targs(:,:)
       integer, allocatable :: itree(:)
       integer *8 ipointer(8), ltree
       integer, allocatable :: ilevel(:)
       integer, allocatable :: nlist1(:),list1(:,:)
       integer, allocatable :: nlist2(:),list2(:,:)
       integer, allocatable :: nlist3(:),list3(:,:)
       integer, allocatable :: nlist4(:),list4(:,:)
       integer mnlist1,mnlist2,mnlist3,mnlist4
       integer, allocatable :: isrcper(:),isrcse(:,:)
       integer, allocatable :: itargper(:),itargse(:,:)
       integer, allocatable :: nlistsrc(:)
       integer, allocatable :: col_ind2(:),row_ind2(:),col_ptr(:)

       allocate(rstmp(ns))

C$OMP PARALLEL DO DEFAULT(SHARED)
       do i=1,ns
         rstmp(i) = 2*rads(i)
       enddo
C$OMP END PARALLEL DO       
       allocate(targs(3,nt))

C$OMP PARALLEL DO DEFAULT(SHARED)       
       do i=1,nt
         targs(1,i) = targets(1,i)
         targs(2,i) = targets(2,i)
         targs(3,i) = targets(3,i)
       enddo
C$OMP END PARALLEL DO       

       isep = 1
       nlmax = 51
       nbmax = 0

       nlevels = 0
       nboxes = 0
       mnlist1 = 0
       mnlist2 = 0
       mnlist3 = 0
       mnlist4 = 0
       mnbors = 27
       ifunif = 0
       nlmin = 0
       iper = 0

       call pts_tree_mem(xyzs,ns,targs,nt,idivflag,ndiv,nlmin,
     1  nlmax,iper,ifunif,nlevels,nboxes,ltree) 

       allocate(centers(3,nboxes),itree(ltree),boxsize(0:nlevels))

       call pts_tree_build(xyzs,ns,targs,nt,idivflag,ndiv,nlmin,
     1   nlmax,iper,ifunif,nlevels,nboxes,ltree,itree,ipointer,centers,
     2   boxsize)
       
       allocate(isrcse(2,nboxes),itargse(2,nboxes))
       allocate(isrcper(ns),itargper(nt))

       call pts_tree_sort(ns,xyzs,itree,ltree,nboxes,nlevels,
     1   ipointer,centers,isrcper,isrcse)
      
       call pts_tree_sort(nt,targs,itree,ltree,nboxes,nlevels,
     1   ipointer,centers,itargper,itargse)

c
c   initialize various tree lists
c
      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0
      mnbors = 27

      isep = 1
      
      call computemnlists(nlevels,nboxes,itree(ipointer(1)),boxsize,
     1  centers,itree(ipointer(3)),itree(ipointer(4)),
     2  itree(ipointer(5)),isep,itree(ipointer(6)),mnbors,
     2  itree(ipointer(7)),iper,mnlist1,mnlist2,mnlist3,mnlist4)
      
      allocate(list1(mnlist1,nboxes),nlist1(nboxes))
      allocate(list2(mnlist2,nboxes),nlist2(nboxes))
      allocate(list3(mnlist3,nboxes),nlist3(nboxes))
      allocate(list4(mnlist4,nboxes),nlist4(nboxes))

      call computelists(nlevels,nboxes,itree(ipointer(1)),boxsize,
     1  centers,itree(ipointer(3)),itree(ipointer(4)),
     2  itree(ipointer(5)),isep,itree(ipointer(6)),mnbors,
     3  itree(ipointer(7)),iper,nlist1,mnlist1,list1,nlist2,
     4  mnlist2,list2,nlist3,mnlist3,list3,
     4  nlist4,mnlist4,list4)

c
c   estimate list of tethered sources
c
c
      allocate(nlistsrc(ns),nlst(nt))
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,ns
        nlistsrc(i) = 0
      enddo
C$OMP END PARALLEL DO      

      nnz = 0
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,nchild,ilev,iss)
C$OMP$PRIVATE(isrc,rtest,ilevup,jbox,i,ncoll,l,lbox,xdis,ydis,zdis)
C$OMP$PRIVATE(llev,distest,itt,itarg,rr,jlev) REDUCTION(+:nnz)
C$OMP$SCHEDULE(DYNAMIC)
      do ibox = 1,nboxes
        nchild = itree(ipointer(4)+ibox-1)
        ilev = itree(ipointer(2)+ibox-1)
        if(nchild.eq.0) then
          do iss=isrcse(1,ibox),isrcse(2,ibox)
             isrc = isrcper(iss)
             rtest = rads(isrc)**2
             ilevup = ceiling(log(rstmp(isrc)/boxsize(ilev))/log(2.0d0))
             jbox = ibox
             ilevup = min(ilevup,ilev)

             do i=1,ilevup
               jbox = itree(ipointer(3)+jbox-1)
             enddo
c
c   loop over colleagues
c
             ncoll = itree(ipointer(6)+jbox-1)
             do l=1,ncoll
               lbox = itree(ipointer(7)+mnbors*(jbox-1)+l-1)
c
c  check if lbox does not intersect with region of interest
c 
               xdis = abs(centers(1,lbox)-xyzs(1,isrc))
               ydis = abs(centers(2,lbox)-xyzs(2,isrc))
               zdis = abs(centers(3,lbox)-xyzs(3,isrc))

               llev = itree(ipointer(2)+lbox-1)
               distest = rads(isrc)+boxsize(llev)/2
               if((xdis.le.distest).and.(ydis.le.distest).and.
     1             (zdis.le.distest)) then
                 
                 do itt=itargse(1,lbox),itargse(2,lbox)
                   itarg = itargper(itt)
                   rr = (xyzs(1,isrc)-targs(1,itarg))**2 + 
     1                (xyzs(2,isrc)-targs(2,itarg))**2 + 
     2                (xyzs(3,isrc)-targs(3,itarg))**2
                   if(rr.le.rtest) then
                     nlistsrc(isrc) = nlistsrc(isrc)+1
                     nnz = nnz + 1
                   endif
                 enddo
               endif
             enddo

c
c   loop over list1 boxes which are larger
c
             jlev = itree(ipointer(2)+jbox-1)
             do l=1,nlist1(jbox)
               lbox = list1(l,jbox)
               llev = itree(ipointer(2)+lbox-1)
               if(llev.lt.jlev) then
                 do itt=itargse(1,lbox),itargse(2,lbox)
                   itarg = itargper(itt)
                   rr = (xyzs(1,isrc)-targs(1,itarg))**2 + 
     1                (xyzs(2,isrc)-targs(2,itarg))**2 + 
     2                (xyzs(3,isrc)-targs(3,itarg))**2
                   if(rr.le.rtest) then
                     nlistsrc(isrc) = nlistsrc(isrc)+1
                     nnz = nnz + 1
                   endif
                 enddo
               endif
             enddo
          enddo
        endif
      enddo
C$OMP END PARALLEL DO      

      allocate(col_ind2(nnz),row_ind2(nnz))
      allocate(col_ptr(ns+1))

      call cumsum(ns,nlistsrc,col_ptr(2))
      col_ptr(1) = 1
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=2,ns+1
        col_ptr(i) = col_ptr(i) + 1
      enddo
C$OMP END PARALLEL DO     

      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,nchild,ilev,iss)
C$OMP$PRIVATE(isrc,rtest,ilevup,jbox,i,ncoll,l,lbox,xdis,ydis,zdis)
C$OMP$PRIVATE(llev,distest,itt,itarg,rr,jlev,nsloc) 
C$OMP$SCHEDULE(DYNAMIC)
      do ibox = 1,nboxes
        nchild = itree(ipointer(4)+ibox-1)
        ilev = itree(ipointer(2)+ibox-1)
        if(nchild.eq.0) then
          do iss=isrcse(1,ibox),isrcse(2,ibox)
             isrc = isrcper(iss)
             nsloc = 0
             rtest = rads(isrc)**2
             ilevup = ceiling(log(rstmp(isrc)/boxsize(ilev))/log(2.0d0))
             jbox = ibox
             ilevup = min(ilevup,ilev)

             do i=1,ilevup
               jbox = itree(ipointer(3)+jbox-1)
             enddo
c
c   loop over colleagues
c
             ncoll = itree(ipointer(6)+jbox-1)
             do l=1,ncoll
               lbox = itree(ipointer(7)+mnbors*(jbox-1)+l-1)
c
c  check if lbox does not intersect with region of interest
c 
               xdis = abs(centers(1,lbox)-xyzs(1,isrc))
               ydis = abs(centers(2,lbox)-xyzs(2,isrc))
               zdis = abs(centers(3,lbox)-xyzs(3,isrc))

               llev = itree(ipointer(2)+lbox-1)
               distest = rads(isrc)+boxsize(llev)/2
               if((xdis.le.distest).and.(ydis.le.distest).and.
     1             (zdis.le.distest)) then
                 
                 do itt=itargse(1,lbox),itargse(2,lbox)
                   itarg = itargper(itt)
                   rr = (xyzs(1,isrc)-targs(1,itarg))**2 + 
     1                (xyzs(2,isrc)-targs(2,itarg))**2 + 
     2                (xyzs(3,isrc)-targs(3,itarg))**2
                   if(rr.le.rtest) then
                     col_ind2(col_ptr(isrc)+nsloc) = isrc
                     row_ind2(col_ptr(isrc)+nsloc) = itarg
                     nsloc = nsloc + 1
                   endif
                 enddo
               endif
             enddo

c
c   loop over list1 boxes which are larger
c
             jlev = itree(ipointer(2)+jbox-1)
             do l=1,nlist1(jbox)
               lbox = list1(l,jbox)
               llev = itree(ipointer(2)+lbox-1)
               if(llev.lt.jlev) then
                 do itt=itargse(1,lbox),itargse(2,lbox)
                   itarg = itargper(itt)
                   rr = (xyzs(1,isrc)-targs(1,itarg))**2 + 
     1                (xyzs(2,isrc)-targs(2,itarg))**2 + 
     2                (xyzs(3,isrc)-targs(3,itarg))**2
                   if(rr.le.rtest) then
                     col_ind2(col_ptr(isrc)+nsloc) = isrc
                     row_ind2(col_ptr(isrc)+nsloc) = itarg
                     nsloc = nsloc + 1
                   endif
                 enddo
               endif
             enddo
          enddo
        endif
      enddo
C$OMP END PARALLEL DO     
      

      if(ns.lt.1000.or.nnz.lt.10000) then
        call conv_to_csc(nnz,nt,col_ind2,row_ind2,row_ptr,col_ind)
      else
        call conv_to_csc_para(nnz,nt,col_ind2,row_ind2,row_ptr,col_ind)
      endif


      return
      end
