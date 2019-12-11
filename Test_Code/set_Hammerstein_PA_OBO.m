function [GAIN PAGLIN] = set_Hammerstein_PA_OBO(input, OBO_FROM_P1DB, pa_coefficients, GINIT, linear_distortion_coefficients)

Q = length(linear_distortion_coefficients);
half_length_Q = (Q + 1) / 2;
N = length(input);
Z0 = 50; %not necessary since calculations are relative but in case you want to know actual power
GAIN = GINIT;

%FIND SS GAIN
lin_input = input * power(10, -30/10); %Take input down by 30 dB to linear point
lin_input_wld = cconv(lin_input, linear_distortion_coefficients);
lin_input_wld = lin_input_wld(1+(half_length_Q):end-(half_length_Q));
ssvout = Memoryless_Polynomial_Amplifier(lin_input_wld, pa_coefficients);
PAOGLIN = 10*log10((ssvout*ssvout')/(N*Z0*0.001));
PAIGLIN = 10*log10((lin_input*lin_input')/(N*Z0*0.001));
PAGLIN = PAOGLIN - PAIGLIN;
%FOUND SS GAIN

%FIND P1DB
temp_LIN = input*GAIN;
temp_LIN_wld = cconv(temp_LIN, linear_distortion_coefficients);
temp_LIN_wld = temp_LIN_wld(1+(half_length_Q):end-(half_length_Q));
temp_NLIN = Memoryless_Polynomial_Amplifier(temp_LIN_wld, pa_coefficients);
TX_POWER_OUT_dBm_LIN = 10*log10((temp_LIN*temp_LIN')/(N*Z0*0.001));
TX_POWER_OUT_dBm_NLIN = 10*log10((temp_NLIN*temp_NLIN')/(N*Z0*0.001));
while TX_POWER_OUT_dBm_NLIN > (TX_POWER_OUT_dBm_LIN + PAGLIN - 1)
   GAIN = GAIN + 0.01;
   temp_LIN = input*GAIN;
   temp_LIN_wld = cconv(temp_LIN, linear_distortion_coefficients);
   temp_LIN_wld = temp_LIN_wld(1+(half_length_Q):end-(half_length_Q));
   temp_NLIN = Memoryless_Polynomial_Amplifier(temp_LIN_wld, pa_coefficients);
   TX_POWER_OUT_dBm_LIN = 10*log10((temp_LIN*temp_LIN')/(N*Z0*0.001));
   TX_POWER_OUT_dBm_NLIN = 10*log10((temp_NLIN*temp_NLIN')/(N*Z0*0.001));
end
P1DB = TX_POWER_OUT_dBm_NLIN;
%FOUND P1DB

if OBO_FROM_P1DB
    %SET OBO GAIN
    while TX_POWER_OUT_dBm_NLIN > (P1DB - OBO_FROM_P1DB)
    GAIN = GAIN - 0.01;
    temp_LIN_wld = cconv(input*GAIN, linear_distortion_coefficients);
    temp_LIN_wld = temp_LIN_wld(1+(half_length_Q):end-(half_length_Q));
    temp_NLIN = Memoryless_Polynomial_Amplifier(temp_LIN_wld, pa_coefficients);
    TX_POWER_OUT_dBm_NLIN = 10*log10((temp_NLIN*temp_NLIN')/(N*Z0*0.001));
    end
end

end