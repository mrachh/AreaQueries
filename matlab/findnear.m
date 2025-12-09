function [rsc] = findnear(srcinfo, rads, targinfo, varargin)
% FINDNEAR    Identify targets close to sources with extent in
%             two and three dimensions.
%
% [rsc] = findnear(srcinfo, rads, targinfo)
% [rsc] = findnear(srcinfo, rads, targinfo, opts)
%
%  Find targets which are contained in |src_{i} - targ_{j}|<= rads_{i}
% 
%  Args:
%  -  srcinfo: double(ndim, ns)  
%       source locations
%  -  rads: double or double(ns)
%       extent of the sources
%  -  targinfo: double(ndimt, nt)
%       target locations
%  -  opts: options struct
%       opts.ifslow (0) use direct algorithm
%       opts.nnz (0) maximum number of source target pairs
%                    expected to be close. If non-zero
%                    then memory estimation routine not 
%                    called
%       opts.ndiv (2) max number of points in a box 
%       opts.idivflag(1) subdivision criterion
%                     0 -> subdivide on sources only
%                     1 -> subdivide on targets only
%                     2 -> subdivide on sources+targets 
%   
%  Returns:
%    * rsc: row sparse compressed struct
%       rsc.row_ptr: (nt+1,1)
%         pointer to col_ind array where list of relevant source
%         patches for target i starts
%       rsc.col_ind: (nnz,1)
%         list of sources relevant for all targets,
%         sorted by target number

  if nargin < 4
    opts = [];
  end

  [ndim, ns] = size(srcinfo);
  if ndim ~= 2 && ndim ~=3 
    error('FINDNEAR: Incorrect leading dimension of srcinfo\n');
  end

  [ndimt, nt] = size(targinfo);
  if ndimt < ndim
    error('FINDNEAR: Leading dimension of targinfo must be same as leading dimension of srcinfo\n');
  end

  tuse = targinfo(1:ndim, :);

  ruse = rads;
  if isscalar(rads)
    ruse = ones(ns,1)*rads;
  end

  ifslow = 0;
  ifpre = 1;
  nnz = 0;
  ndiv = 2;
  idivflag = 0;
  if isfield(opts, 'ifslow'), ifslow = opts.ifslow; end
  if isfield(opts, 'ndiv'), ndiv = opts.ndiv; end
  if isfield(opts, 'idivflag'), idivflag = opts.idivflag; end
  if isfield(opts, 'nnz')
    ifpre = 0;
    nnz = opts.nnz;
  end

  rsc = [];
  ntp1 = nt+1;
  row_ptr = zeros(ntp1,1);
  if ndim == 2
    if ifpre
      if ifslow
        mex_id_ = 'findnearslowmem2d(c i double[xx], c i int[x], c i double[x], c i double[xx], c i int[x], c io int[x])';
[nnz] = areaqueries(mex_id_, srcinfo, ns, ruse, tuse, nt, nnz, 2, ns, 1, ns, 2, nt, 1, 1);
      else
        mex_id_ = 'findnearmem2d(c i double[xx], c i int[x], c i double[x], c i double[xx], c i int[x], c io int[x])';
[nnz] = areaqueries(mex_id_, srcinfo, ns, ruse, tuse, nt, nnz, 2, ns, 1, ns, 2, nt, 1, 1);
      end
    end
    col_ind = zeros(nnz,1);
    if ifslow
      mex_id_ = 'findnearslow2d(c i double[xx], c i int[x], c i double[x], c i double[xx], c i int[x], c io int[x], c io int[x])';
[row_ptr, col_ind] = areaqueries(mex_id_, srcinfo, ns, ruse, tuse, nt, row_ptr, col_ind, 2, ns, 1, ns, 2, nt, 1, ntp1, nnz);
    else
      mex_id_ = 'findnear2d(c i double[xx], c i int[x], c i double[x], c i double[xx], c i int[x], c io int[x], c io int[x])';
[row_ptr, col_ind] = areaqueries(mex_id_, srcinfo, ns, ruse, tuse, nt, row_ptr, col_ind, 2, ns, 1, ns, 2, nt, 1, ntp1, nnz);
    end

  else
    if ifpre
      if ifslow
        mex_id_ = 'findnearslowmem3d(c i double[xx], c i int[x], c i double[x], c i int[x], c i double[xx], c i int[x], c io int[x])';
[nnz] = areaqueries(mex_id_, srcinfo, ns, ruse, ndim, tuse, nt, nnz, 3, ns, 1, ns, 1, 3, nt, 1, 1);
      else
        mex_id_ = 'findnearmem3d(c i double[xx], c i int[x], c i double[x], c i int[x], c i double[xx], c i int[x], c io int[x])';
[nnz] = areaqueries(mex_id_, srcinfo, ns, ruse, ndim, tuse, nt, nnz, 3, ns, 1, ns, 1, 3, nt, 1, 1);
      end
    end
    col_ind = zeros(nnz,1);
    if ifslow
      mex_id_ = 'findnearslow3d(c i double[xx], c i int[x], c i double[x], c i int[x], c i double[xx], c i int[x], c io int[x], c io int[x])';
[row_ptr, col_ind] = areaqueries(mex_id_, srcinfo, ns, ruse, ndim, tuse, nt, row_ptr, col_ind, 3, ns, 1, ns, 1, 3, nt, 1, ntp1, nnz);
    else
      mex_id_ = 'findnear3d(c i double[xx], c i int[x], c i double[x], c i int[x], c i double[xx], c i int[x], c io int[x], c io int[x])';
[row_ptr, col_ind] = areaqueries(mex_id_, srcinfo, ns, ruse, ndim, tuse, nt, row_ptr, col_ind, 3, ns, 1, ns, 1, 3, nt, 1, ntp1, nnz);
    end
  end
  rsc.row_ptr = row_ptr;
  rsc.col_ind = col_ind;

end
