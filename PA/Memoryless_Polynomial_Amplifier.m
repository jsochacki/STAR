function [y] = Memoryless_Polynomial_Amplifier(x, coefficients)
% The equation is
% y1=c((N+1)/2)*x1*abs(x1)^(N-1) + c(((N+1)/2) - 2)*x1*abs(x1)^(N-3) + ... + c1*x1
% where
% yn=c((N+1)/2)*xn*abs(xn)^(N-1) + c(((N+1)/2) - 2)*xn*abs(xn)^(N-3) + ... + c1*xn
%
% In matrix form Y = X*coefficients
% | y1 | _ |  x1*(abs(x1)^(N-1)) x1*(abs(x1)^(N-3)) ... x1  | | c((N+1)/2)       |
% | y2 | - |  x2*(abs(x2)^(N-1)) x2*(abs(x2)^(N-3)) ... x2  | | c(((N+1)/2) - 2) |
%  .
%  .
% | yn | - |  xn*(abs(xn)^(N-1)) xn*(abs(xn)^(N-3)) ... xn  | | c1               |
%
% Y : These are the output points for the local interpolation
% x : %These are the input points for the local interpolation
% N : %This is the order of the polynomial.  It should be equal to the
%     highest odd nonlinear term that you need to generate a value for

% Turn inputs so we can do matrix generation and inversion
trn1=0;
if size(x,2)> size(x,1), trn1=1; x=x.';, end;
if size(coefficients,2)> size(coefficients,1), coefficients=coefficients.';, end;

N = length(coefficients) - 1;

% Generate x matrix and pick odd coefficients only
X = zeros(length(x), (N + 1) / 2);
odd_coefficients = zeros((N + 1) / 2, 1);
column = 1;
for n = 1:2:N
   X(:, column) = x.*(abs(x).^(N - n));
   odd_coefficients(column, 1) = coefficients(n);
   column = column + 1;
end

% Least Squares Inversion
y = X * odd_coefficients;

if trn1, y=y.';, end;
end