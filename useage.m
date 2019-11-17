clear all
clc

x = [1 2 3 4 5 6 7 8 9 10];
y = [2 4 6 7.5 9 10 11 11.3 11.7 12] + j*[0.001 0.01 0.03 0.06 0.1 0.2 0.3 0.4 0.4 0.3];
coefficients = Least_Squares_Odd_Polynomial_Solution(y, x, 7);
y2 = Memoryless_Polynomial_Amplifier(x, coefficients);
GAIN = 2;

%AM - AM replication
figure(1)
plot(10*log10(x.*x'.'), 10*log10(y.*y'.'))
hold on
plot(10*log10(x.*x'.'), 10*log10(y2.*y2'.'))

figure(2)
%AM - PM replication
plot(10*log10(x.*x'.'), atan(imag(y)./real(y)))
hold on
plot(10*log10(x.*x'.'), atan(imag(y2)./real(y2)))


extended_x = 0:0.01:20;
y3 = Memoryless_Polynomial_Amplifier(extended_x, coefficients);

%AM - AM replication
figure(3)
plot(10*log10(x.*x'.'), 10*log10(y.*y'.'))
hold on
plot(10*log10(extended_x.*extended_x'.'), 10*log10(y3.*y3'.'))
plot(10*log10(extended_x.*extended_x'.'), 10*log10((GAIN*extended_x).*(GAIN*extended_x)'.'))

figure(4)
%AM - PM replication
plot(10*log10(x.*x'.'), atan(imag(y)./real(y)))
hold on
plot(10*log10(extended_x.*extended_x'.'), atan(imag(y3)./real(y3)))

figure(5)
plot(10*log10(x.*x'.'), 10*log10(y2.*y2'.') - 10*log10(x.*x'.'))
hold on
plot(10*log10(extended_x.*extended_x'.'), 10*log10(y3.*y3'.') - 10*log10(extended_x.*extended_x'.'))