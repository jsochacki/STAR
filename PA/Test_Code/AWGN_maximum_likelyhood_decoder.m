function [decoded_complex_stream] = AWGN_maximum_likelyhood_decoder(received_baseband_data, Constellation, complex_mapping)
trn=0;
if size(received_baseband_data,1) > size(received_baseband_data,2), received_baseband_data=received_baseband_data.';, trn=1;, end;

N = length(received_baseband_data);
decoded_complex_stream = zeros(1, N);

for n = 1:1:length(received_baseband_data)
    temp = received_baseband_data(1, n) - complex_mapping;
    EUCDIS = temp.*temp'.';
    [MAG IND] = sort(EUCDIS, 'ascend');
    if sum(MAG(1) == MAG(2:1:length(MAG)))
        decoded_complex_stream(n) = randsrc(1, 1, Constellation(EUCDIS == MAG(1)));
    else
        decoded_complex_stream(n) = Constellation(IND(1));
    end
end

if trn, decoded_complex_stream = decoded_complex_stream.';, end;

end