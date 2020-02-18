function [SNR_dB] = Measure_SNR(rx_symbols, noiseless_symbols)

N = length(rx_symbols);
if size(rx_symbols,1)> size(rx_symbols,2), rx_symbols=rx_symbols.';, end;
if size(noiseless_symbols,1)> size(noiseless_symbols,2), noiseless_symbols=noiseless_symbols.';, end;

VERROR = rx_symbols - noiseless_symbols;
PERROR = (1/N)*(VERROR * VERROR') ;
PREF = (1/N)*(noiseless_symbols * noiseless_symbols') ;
SNR_dB = 10*log10(PREF / PERROR);

end