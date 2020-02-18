function [result] = matrix_convolve_and_truncate(x, b)
% Turn inputs so we can do matrix generation and inversion
trn1=0;
if size(x,2)> size(x,1), trn1=1; x=x.'; b=b.';, end;

N = length(b);
% Generate X matrix
X = zeros(N, (length(x) - N) + 1);
for n = 1:1:N
   X(n, :) = x(n:1:(end-(N-n))).';
end

X = X.';

result = X * flipud(b);

if trn1, result=result.';, end;

end