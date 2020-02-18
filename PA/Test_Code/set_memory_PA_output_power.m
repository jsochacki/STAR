function [GAIN] = set_memory_PA_output_power(input, output_power_dBm, pa_coefficients, GINIT, PA_ORDER, Q)

N = length(input);
Z0 = 50; %not necessary since calculations are relative but in case you want to know actual power
GAIN = GINIT;

P_TARGET = Z0 * 0.001 * power(10, output_power_dBm / 10);
%FIND PTARGET
temp_NLIN = Memory_Polynomial_Amplifier(input*GAIN, pa_coefficients, PA_ORDER, Q);
TX_POWER_OUT_dBm_NLIN = 10*log10((temp_NLIN*temp_NLIN')/(N*Z0*0.001));
while P_TARGET > TX_POWER_OUT_dBm_NLIN
   GAIN = GAIN + 0.01;
   temp_NLIN = Memory_Polynomial_Amplifier(input*GAIN, pa_coefficients, PA_ORDER, Q);
   TX_POWER_OUT_dBm_NLIN = 10*log10((temp_NLIN*temp_NLIN')/(N*Z0*0.001));
end
%FOUND PTARGET

end