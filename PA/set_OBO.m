function [GAIN PAGLIN] = set_OBO(input, OBO_FROM_P1DB, pa_coefficients, GINIT)

N = length(input);
Z0 = 50; %not necessary since calculations are relative but in case you want to know actual power
GAIN = GINIT;

%FIND SS GAIN
lin_input = input * power(10, -30/10); %Take input down by 30 dB to linear point
ssvout = Memoryless_Polynomial_Amplifier(lin_input, pa_coefficients);
PAOGLIN = 10*log10((ssvout*ssvout')/(N*Z0*0.001));
PAIGLIN = 10*log10((lin_input*lin_input')/(N*Z0*0.001));
PAGLIN = PAOGLIN - PAIGLIN;
%FOUND SS GAIN

%FIND P1DB
temp_LIN = input*GAIN;
temp_NLIN = Memoryless_Polynomial_Amplifier(temp_LIN, pa_coefficients);
TX_POWER_OUT_dBm_LIN = 10*log10((temp_LIN*temp_LIN')/(N*Z0*0.001));
TX_POWER_OUT_dBm_NLIN = 10*log10((temp_NLIN*temp_NLIN')/(N*Z0*0.001));
while TX_POWER_OUT_dBm_NLIN > (TX_POWER_OUT_dBm_LIN + PAGLIN - 1)
   GAIN = GAIN + 0.01;
   temp_LIN = input*GAIN;
   temp_NLIN = Memoryless_Polynomial_Amplifier(temp_LIN, pa_coefficients);
   TX_POWER_OUT_dBm_LIN = 10*log10((temp_LIN*temp_LIN')/(N*Z0*0.001));
   TX_POWER_OUT_dBm_NLIN = 10*log10((temp_NLIN*temp_NLIN')/(N*Z0*0.001));
end
P1DB = TX_POWER_OUT_dBm_NLIN;
%FOUND P1DB

if OBO_FROM_P1DB
    %SET OBO GAIN
    while TX_POWER_OUT_dBm_NLIN > (P1DB - OBO_FROM_P1DB)
    GAIN = GAIN - 0.01;
    temp_NLIN = Memoryless_Polynomial_Amplifier(input*GAIN, pa_coefficients);
    TX_POWER_OUT_dBm_NLIN = 10*log10((temp_NLIN*temp_NLIN')/(N*Z0*0.001));
    end
end

end