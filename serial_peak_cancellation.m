function [processed_signal current_PAPR_dB] = serial_peak_cancellation(input_signal, filter_h, target_PAPR_dB, iteration_count)

processed_signal = input_signal;
impulse = filter_h ./ max(abs(filter_h));
target_PAPR = power(10, target_PAPR_dB / 10);
N = length(processed_signal);
impulse_half_length = (length(impulse) - 1) / 2;
[current_PAPR_dB current_PAPR] = PAPR_dB(processed_signal, []);

iteration = 1;
while iteration_count >= iteration
   
   magnitude = sqrt(processed_signal.*processed_signal'.');
   current_peak_magnitude = max(magnitude);
   index = find(current_peak_magnitude == magnitude);
   
   if isempty(index)
      break
   end
   
   current_average_power = (processed_signal*processed_signal') / N;
   current_target_peak_power = current_average_power * target_PAPR;
   current_target_peak_voltage = sqrt(current_target_peak_power);
   
   indicies = ((index - impulse_half_length):1:(index + impulse_half_length));

   amplitude_scaling = current_peak_magnitude - current_target_peak_power;
   phase_scaling = angle(processed_signal(index));
   
   correction_pulse = -(impulse .* amplitude_scaling .* exp(j * phase_scaling));
   
   if (sum(indicies > 0) == (2*impulse_half_length + 1)) ...
       && (sum(indicies <= N) == (2*impulse_half_length + 1))
      processed_signal(indicies) = processed_signal(indicies) + correction_pulse;
   elseif (indicies(1) < 0)
      samples_to_drop = 1 - indicies(1);
      processed_signal(indicies(1:1:(end-samples_to_drop))) = processed_signal(indicies(1:1:(end-samples_to_drop))) + correction_pulse((1 + samples_to_drop):1:end);
   else
      samples_to_drop = indicies(end) - N;
      processed_signal(indicies(1:1:(end-samples_to_drop))) = processed_signal(indicies(1:1:(end-samples_to_drop))) + correction_pulse(1:1:end-samples_to_drop);
   end
   
   iteration = iteration + 1;
end

current_PAPR_dB = PAPR_dB(processed_signal, []);
   
end