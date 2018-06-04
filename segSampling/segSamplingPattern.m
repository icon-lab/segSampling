function [minIntrVec,stat,N] = segSamplingPattern(pdf,iter,tol, initialPDF)
%	A slightly modified version of Michael Lustig's samplingPattern for generating
%	segregated masks
%
% a monte-carlo algorithm to generate a sampling pattern with 
% minimum peak interference. The number of samples will be
% sum(pdf) +- tol
%
%	pdf - probability density function to choose samples from
%	iter - number of tries
%	tol  - the deviation from the desired number of samples in samples
%   initialPDF - initial PDF is used for interference calculations
%
% returns:
%	mask - sampling pattern
%	stat - vector of min interferences measured each try
%	N    - index in which the minimum interference occured

h = waitbar(0);

pdf(find(pdf>1)) = 1;
K = sum(pdf(:));

minIntr = 1e99;
minIntrVec = zeros(size(pdf));
for n=1:iter
	tmp = zeros(size(pdf));
	while abs(sum(tmp(:)) - K) > tol
		tmp = rand(size(pdf))<pdf;
	end
	x = (tmp./initialPDF);
	x (isnan(x(:))) = 0;
	TMP = ifft2(x);
	if max(abs(TMP(2:end))) < minIntr
		minIntr = max(abs(TMP(2:end)));
		minIntrVec = tmp;
		N = n;
	end
	stat(n) = max(abs(TMP(2:end)));
	waitbar(n/iter,h);
end
close(h);


