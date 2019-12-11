function [coefficients_out X] = Hammerstein_Memoryless_Odd_Polynomial_Solution(Y, x, N, linear_coefficients)
% Turn inputs so we can do matrix generation and inversion
trn1=0;
if size(Y,2)> size(Y,1), Y=Y.';, end;
if size(x,2)> size(x,1), trn1=1; x=x.';, end;

% Generate x matrix
X = zeros(length(x), (N + 1) / 2);
column = 1;
for n = 1:2:N
   X(:, column) = x.*(abs(x).^(N - n));
   column = column + 1;
end

%Dont zero pad like you do in the amplifier, here you only have the x you
%have (here the x will typically actually be output values since this is a
%solution) so you need to have a truncated U matrix to solve for all x's you
%actually have so you need to use the truncated convolution
%Get initial memoryless pd coefficients from signal with memory
U = [];
for column = 1:1:((N + 1) / 2)
   U = [U matrix_convolve_and_truncate(X(:, column), linear_coefficients.')];
end

Q = length(linear_coefficients);
Y = Y((1+((Q - 1) / 2)):1:(end-((Q - 1) / 2)));

coefficients = pinv(U)*Y;
%coefficients = lsqminnorm(U, Y);

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

if trn1, X=X.';, end;
if trn1, coefficients_out=coefficients_out.';, end;
end