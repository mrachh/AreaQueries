c     test driver for findnear2d.
c     compares fast to slow (naive) method for a bunch of random disks
c     and targets. The sorted disk indices for each target for the two
c     methods are compared.

      implicit real *8 (a-h,o-z)
      real *8, allocatable :: src(:,:),targ(:,:),rads(:)
      integer, allocatable :: row_ptr(:),row_ptr2(:)
      integer, allocatable :: col_ind(:),col_ind2(:)
      integer, allocatable :: isort(:),wsort(:),isort2(:),wsort2(:)
c     for timing
      integer crate,t1,t2,t3
      real *8 t,ts
      
      call prini(6,13)


      ns = 10001
      nt = 50000
      
      allocate(src(2,ns),targ(2,nt),rads(ns))
      allocate(row_ptr(nt+1),row_ptr2(nt+1))
      allocate(isort(nt),wsort(nt),isort2(nt),wsort2(nt))

c     ball centers in unit square with random radii
      do i=1,ns
        src(1,i) = hkrand(0)
        src(2,i) = hkrand(0)
        rads(i) = 2.0d0**(-10*hkrand(0))/100.0d0
      enddo

      do i=1,nt
        targ(1,i) = hkrand(0)
        targ(2,i) = hkrand(0)
      enddo

      nnz = 0
      nnz2 = 0
      call system_clock(t1)
      call findnearslowmem2d(src,ns,rads,targ,nt,nnz2)
      call system_clock(t2,crate)
      ts = (t2-t1)/float(crate)
      call findnearmem2d(src,ns,rads,targ,nt,nnz)
      call system_clock(t3)
      t = (t3-t2)/float(crate)
      
      print *,'ns=',ns,' nt=',nt,' nnz=',nnz,' nnz/nt=',float(nnz)/nt
      print '("tmem slow:  ",f6.3," sec,  tmem fast:  ",f6.3," sec")',
     1    ts,t

      
      if(nnz.ne.nnz2) then
        call prinf('number of non zero elements dont match*',i,0)
        stop
      endif

      allocate(col_ind(nnz),col_ind2(nnz))
      call system_clock(t1)
      call findnearslow2d(src,ns,rads,targ,nt,row_ptr2,col_ind2)
      call system_clock(t2,crate)
      ts = (t2-t1)/float(crate)
      call findnear2d(src,ns,rads,targ,nt,row_ptr,col_ind)
      call system_clock(t3)
      t = (t3-t2)/float(crate)
      print '("tfind slow: ",f6.3," sec,  tfind fast: ",f6.3," sec")',
     1    ts,t
      
      do i=1,nt
        n1 = row_ptr(i+1)-row_ptr(i)
        n2 = row_ptr2(i+1)-row_ptr2(i)
        if(n1.ne.n2) then
          call prinf('number of sources doesnt match for target i=*',
     1       i,1)
          stop
        endif

        call sorti(n1,col_ind(row_ptr(i)),wsort)
        call sorti(n2,col_ind2(row_ptr2(i)),wsort2)

        erra = 0
        do j=1,n1
          isort(j) = col_ind(row_ptr(i)+wsort(j)-1)
          isort2(j) = col_ind2(row_ptr2(i)+wsort2(j)-1)
          erra = erra + abs(isort(j)-isort2(j))
        enddo

        if(erra.ne.0) then
          call prinf('list of sources dont match for target i=*',i,1)
          call prinf('correct sources=*',isort2,n2)
          call prinf('computed sources=*',isort,n1)
          stop
        endif
      enddo
      
      print *, "2D test successfully completed"



      stop
      end
