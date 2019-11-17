function [coefficients X] = Memory_Polynomial_Solution(y, x, N, Q)
% Turn inputs so we can do matrix generation and inversion
trn1=0;
if size(y,2)> size(y,1), y=y.';, end;
if size(x,2)> size(x,1), trn1=1; x=x.';, end;

%Dont zero pad like you do in the amplifier, here you only have the y you
%have so you need to have a truncated x matrix to solve for all y's you
%actually have
% Generate X matrix
X = zeros(((length(x) - Q) + 1), ((N + 1) / 2) * Q);
SUBMATRIX = zeros(Q, (length(x) - Q) + 1);
term = 1;
for k = N:-2:1
   for n = 1:1:Q
      SUBMATRIX(n, :) = (x(n:1:(end-(Q-n))).').*(abs(x(n:1:(end-(Q-n))).').^(N - k));
   end
   
   X(1:1:((length(x) - Q) + 1), (Q*(term-1)+1):1:(Q*term)) =  SUBMATRIX.';
   term = term + 1;
end

%Generate Y matrix
Y = y((1+((Q - 1) / 2)):1:(end-((Q - 1) / 2)));

% Least Squares Inversion
%Dont try and use inv for linear systems in matlab, it is usually singular
%due to the LU decomposition
%coefficients = inv(X'*X) * X' * Y;
%In addition in the Rank of X is less than the number of columns of A (i.e.
%non square matrix) then A\B is not the minimum norm solution and you
%can compute the minimum norm least-squares solution using
%a = lsqminnorm(X,Y) (COD complete orthogonal decomposition) or
%a = pinv(X)*Y (SVD)
%coefficients = X\Y;
coefficients = pinv(X)*Y;
%coefficients = lsqminnorm(X, Y);


if trn1, X=X.';, end;
if trn1, coefficients=coefficients.';, end;
end