ns = 10001;
nt = 50001;
srcinfo = rand(2, ns);
targinfo = rand(2, nt);
rads = 0.1*rand(ns,1);

tic, rsc = findnear(srcinfo, rads, targinfo); toc;

% 3d test
srcinfo = rand(3, ns);
targinfo = rand(3, nt);
rads = 0.1*rand(ns,1);

tic, rsc = findnear(srcinfo, rads, targinfo); toc;
