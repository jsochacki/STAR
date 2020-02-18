function [y] = Memory_Polynomial_Amplifier(x, coefficients, N, Q)

% Turn inputs so we can do matrix generation and inversion
trn1=0;
if size(x,2)> size(x,1), trn1=1; x=x.';, end;
if size(coefficients,2)> size(coefficients,1), coefficients=coefficients.';, end;

% Generate X matrix
X = zeros(((length(x) + Q) - 1), ((N + 1) / 2) * Q);
x_padded = [zeros(Q - 1, 1); x; zeros(Q - 1, 1)];
SUBMATRIX = zeros(Q, (length(x) + Q) - 1);
term = 1;
for k = N:-2:1
   for n = 1:1:Q
      SUBMATRIX(n, :) = (x_padded(n:1:(end-(Q-n))).').*(abs(x_padded(n:1:(end-(Q-n))).').^(N - k));
   end
   
   X(1:1:((length(x) + Q) - 1), (Q*(term-1)+1):1:(Q*term)) =  SUBMATRIX.';
   term = term + 1;
end

% Least Squares Inversion
y = X * coefficients;

y = y((1+((Q - 1) / 2)):end-((Q - 1) / 2));

if trn1, y=y.';, end;
end