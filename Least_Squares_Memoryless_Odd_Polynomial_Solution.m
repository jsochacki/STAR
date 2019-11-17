function [coefficients_out X] = Least_Squares_Memoryless_Odd_Polynomial_Solution(Y, x, N)
% The equation is
% y1=c((N+1)/2)*x1^N + c(((N+1)/2) - 1)*x1^(N-2) + ... + c1*x1
% where
% yn=c((N+1)/2)*xn^N + c(((N+1)/2) - 1)*xn^(N-2) + ... + c1*xn
%
% In matrix form Y = X*coefficients
% | y1 | _ |  x1^N x1^(N-2) ... x1  | | c((N+1)/2)       |
% | y2 | - |  x2^N x2^(N-2) ... x2  | | c(((N+1)/2) - 1) |
%  .
%  .
% | yn | - |  xn^N xn^(N-2) ... xn  | | c1               |
%
% Y : These are the output points for the local interpolation
% x : %These are the input points for the local interpolation
% N : %This is the order of the polynomial.  It should be equal to the
%     highest odd nonlinear term that you need to generate a value for

% Turn inputs so we can do matrix generation and inversion
trn1=0; trn2=0;
if size(Y,2)> size(Y,1), trn1=1; Y=Y.';, end;
if size(x,2)> size(x,1), trn2=1; x=x.';, end;

% Generate x matrix
X = zeros(length(x), (N + 1) / 2);
column = 1;
for n = 1:2:N
   X(:, column) = x.*(abs(x).^(N - n));
   column = column + 1;
end

% Least Squares Inversion
%coefficients = inv(X'*X) * X' * Y;
coefficients = pinv(X)*Y;
%coefficients = lsqminnorm(X, Y);

%Make coefficients vector full polynomial
coefficients_out = zeros(N, 1);
index = 1;
for n = 1:1:(N+1)
   if mod(n, 2)
      coefficients_out(n) = coefficients(index);
      index = index + 1;
   else
      coefficients_out(n) = 0;
   end
end

if trn2, X=X.';, end;
if trn2, coefficients_out=coefficients_out.';, end;
end