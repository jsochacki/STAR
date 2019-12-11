function [output_signal, sigma] = generate_awgn_from_EsNo(input_signal, sigma_in, EsNo_dB, oversampling_rate)
%generate_awgn Makes noise signal based on EsNo input in dB

   N = length(input_signal);
   if (nargin == 2)
       sigma = sigma_in;
   elseif (nargin == 4)
       EsNo = power(10, EsNo_dB / 10);
       mup2 = (1/N) * (input_signal*input_signal');
       mup2 = mup2 * oversampling_rate;
       sigma = sqrt(mup2 / (2*EsNo));
   end

   output_signal = sigma*randn(1, N) + j*sigma*randn(1, N);
end
