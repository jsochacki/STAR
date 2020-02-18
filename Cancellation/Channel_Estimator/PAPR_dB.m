function [PAPR_dB PAPR CF] = PAPR_dB(waveform, symbols_and_pdfs)
%%%symbols_and_pdfs is a two column vector where the rows are the symbols
%%%and their corresponding pdf.  The first column should be the symbols
%%%the second column should be the pdf
%%%%i.e. [0.707+j*0.707 52;0.707-j*0.707 55;-0.707+j*0.707 12;-0.707-j*0.707 32]
%%%if there are no stastics just use [] for symbols_and_pdfs
if ~isempty(symbols_and_pdfs)
    symbols_and_probabilities = []; symbols_and_probabilities(:,1) = symbols_and_pdfs(:,1);
    for n = 1:1:length(symbols_and_pdfs(:,2))
        symbols_and_probabilities(n,2) = symbols_and_pdfs(n,2) / sum(symbols_and_pdfs(:,2));
    end
    powers_and_probabilities = [symbols_and_probabilities(:,1).*symbols_and_probabilities(:,1)'.' symbols_and_probabilities(:,2)];
    PAPR = max(waveform.*waveform'.') / (sum(prod(powers_and_probabilities,2)));
 else
    PAPR = max(waveform.*waveform'.') / (sum(waveform.*waveform'.') / length(waveform));
end

PAPR_dB = 10*log10(PAPR);
CF = sqrt(power(10, PAPR/10));

end