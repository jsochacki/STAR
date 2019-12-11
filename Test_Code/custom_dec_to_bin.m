function [binary]=custom_dec_to_bin(decimal_stream,BITS_PER_WORD)
    trn=0;
    binary=[]; binary_table=power(2,(BITS_PER_WORD-1):-1:0);
    if size(decimal_stream,1) > 1, trn=1; decimal_stream=decimal_stream.';, end;
    for nn=1:1:length(decimal_stream)
        Symbol=decimal_stream(nn); bits=[];
        for n=1:1:BITS_PER_WORD
            bits(n)=(Symbol/(binary_table(n))) >= 1;
            Symbol=Symbol-(binary_table(n))*bits(n);
        end
        binary=[binary bits];
    end
    if trn, binary=binary.';, end;
end