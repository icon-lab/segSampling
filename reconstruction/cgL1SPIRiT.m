function [res, RESVEC] = cgL1SPIRiT(DATA, GOP, nIter, lambda, wavWeight, nIterSplit)

if nargin < 5
        wavWeight = 0;
        nIterSplit = 1;
end

[sx, sy, N] = size(DATA);
ssx = 2^ceil(log2(sx)); 
ssy = 2^ceil(log2(sy));
ss = max(ssx, ssy);
W = Wavelet('Daubechies',4,4);

kernel = getKernel(GOP);
kSize = [size(kernel,1),size(kernel,2)];

idx_acq = find(abs(DATA)>0);
idx_nacq = find(abs(DATA)==0);
empNum = length(idx_nacq(:));

yy = GOP*DATA; yy = [-yy(:); idx_nacq(:)*0];

xn = zeros(empNum,1);
for i=1:nIterSplit
    [xn,FLAG,RELRES,ITER,RESVEC] = lsqr(@aprod,yy,1e-6,nIter, speye(empNum,empNum),speye(empNum,empNum),xn,GOP,sx,sy,N,idx_nacq, lambda);
    if wavWeight > 0
        X = DATA;
        X(idx_nacq) = xn;
        X = ifft2c(X);
        X = zpad(X,ss,ss,N);
        X = reshape(X,[ss,ss,N]);
        X = W*(X);
        X = SoftThresh(X, wavWeight);
        X = W'*(X);
        X = reshape(X,[ss,ss,N]);
        X = crop(X,sx,sy,N);
        X = fft2c(X);
        xn = X(idx_nacq);
    end
end
res = DATA;
res(idx_nacq) = xn;


function [res,tflag] = aprod(x,GOP,sx,sy,nCoils,idx_nacq, lambda,tflag)

	kernel = getKernel(GOP);
	kSize = [size(kernel,1),size(kernel,2)];

	if strcmp(tflag,'transp');
		tmpy = reshape(x(1:sx*sy*nCoils),sx,sy,nCoils);
        res = GOP'*tmpy;
        res = res(idx_nacq)+ x(sx*sy*nCoils+1:end)*lambda;
	
    else
		tmpx = zeros(sx,sy,nCoils);
		tmpx(idx_nacq) = x;
		res = GOP*tmpx;
		res = [res(:) ; lambda*x(:)];
    end


