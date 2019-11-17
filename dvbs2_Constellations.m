function [Complex_Alphabet Binary_Alphabet Decimal_Alphabet BITS_PER_WORD] = dvbs2_Constellations(MODCOD)
%%%% The following constellations are designed and ordered according to the DVB-S2 standard, ETSI EN 302 307 V1.2.1 (2009-08)
%%%% THEY ARE IN ORDER SO alphabet(1) is the complex point
%%%% associated with a binary 1, alphabet(2) is the complex point
%%%% associated with a binary 2 ect .. for QPSK, 8PSK, and 16APSK

if MODCOD == 1
    BITS_PER_WORD = 2;
    CONSTELLATION_NAME = 'QPSK';
    Complex_Alphabet = exp(1j*([0,3,1,2].'*pi/2+pi/4)).';
    Decimal_Alphabet = [0 1 2 3];
    Binary_Alphabet = custom_dec_to_bin(Decimal_Alphabet,BITS_PER_WORD);
elseif MODCOD == 2
    BITS_PER_WORD = 3;
    CONSTELLATION_NAME = '8PSK';
    Complex_Alphabet = exp(1j*([0,7,3,4,1,2,6,5].'*pi/4+pi/4)).';
    Decimal_Alphabet = [0 1 2 3 4 5 6 7];
    Binary_Alphabet = custom_dec_to_bin(Decimal_Alphabet,BITS_PER_WORD);
elseif (MODCOD >= 3 && MODCOD <= 8) || MODCOD == 18
    CONSTELLATION_NAME = '16APSK';
    if MODCOD == 3, RingRatio = 3.15; %outer/inner
    elseif MODCOD == 4, RingRatio = 2.85;
    elseif MODCOD == 5, RingRatio = 2.75;
    elseif MODCOD == 6, RingRatio = 2.70;
    elseif MODCOD == 7, RingRatio = 2.60;
    elseif MODCOD == 8, RingRatio = 2.57;
    elseif MODCOD == 18, RingRatio = 3.15;
    end
    BITS_PER_WORD = 4;
    Complex_Alphabet = [RingRatio*exp(1j*([0,9,3,6,11,10,4,5,1,8,2,7].'*pi/6+pi/4)); exp(1j*([0,3,1,2].'*pi/2+pi/4))].';
    Decimal_Alphabet = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
    Binary_Alphabet = custom_dec_to_bin(Decimal_Alphabet,BITS_PER_WORD);
elseif MODCOD >= 9 && MODCOD <= 15
    CONSTELLATION_NAME = '32APSK';
    if MODCOD == 9, RingRatio(1) = 2.84; RingRatio(2) = 5.27; %middle/inner & outer/inner
    elseif MODCOD == 10, RingRatio(1) = 2.72; RingRatio(2) = 4.87;
    elseif MODCOD == 11, RingRatio(1) = 2.64; RingRatio(2) = 4.64;
    elseif MODCOD == 12, RingRatio(1) = 2.54; RingRatio(2) = 4.33;
    elseif MODCOD == 13, RingRatio(1) = 2.53; RingRatio(2) = 4.30;
    elseif MODCOD == 14, RingRatio(1) = 1.95; RingRatio(2) = 3.03; %trial modcod (remove, not part of standard)
    elseif MODCOD == 15, RingRatio(1) = 2.84*1; RingRatio(2) = 5.27*0.8; %trial modcod (remove, not part of standard)
    end
    BITS_PER_WORD = 5;
    Complex_Alphabet = [RingRatio(1)*exp(1j*([0,1,9,8,3,2,6,7].'*pi/6+pi/4)); RingRatio(2)*exp(1j*([0,2,13,11,5,3,8,10].'*pi/8+pi/8)); %pt0-15
        RingRatio(1)*exp(1j*(11*pi/6+pi/4)); exp(1j*(0*pi/2+pi/4)); RingRatio(1)*exp(1j*(10*pi/6+pi/4)); exp(1j*(3*pi/2+pi/4)); %pt16-19
        RingRatio(1)*exp(1j*(4*pi/6+pi/4)); exp(1j*(1*pi/2+pi/4)); RingRatio(1)*exp(1j*(5*pi/6+pi/4)); exp(1j*(2*pi/2+pi/4)); %pt20-23
        RingRatio(2)*exp(1j*([15,1,14,12,6,4,7,9].'*pi/8+pi/8))].'; %pt24-31
    Decimal_Alphabet = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31];
elseif MODCOD == 16 %%NOT PART OF DVBS2.  16 QAM not designed to ANY standard.
    BITS_PER_WORD = 4;
    CONSTELLATION_NAME = '16_QAM';
    Complex_Alphabet = [];
    for n = 0:1:(nextpow2(2^BITS_PER_WORD)-1)
        Complex_Alphabet = [Complex_Alphabet exp(j*pi/4).*ones(1,nextpow2(2^BITS_PER_WORD))-j*n*(2*real(exp(j*pi/4))/(nextpow2(2^BITS_PER_WORD)-1))-[(cumsum(ones(1,nextpow2(2^BITS_PER_WORD)))-1).*(2*real(exp(j*pi/4))/(nextpow2(2^BITS_PER_WORD)-1))]];
    end
    Decimal_Alphabet = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
    Binary_Alphabet = custom_dec_to_bin(Decimal_Alphabet,BITS_PER_WORD);
elseif MODCOD == 17 %%NOT PART OF DVBS2.  32 QAM not designed to ANY standard.
    BITS_PER_WORD = 5;
    CONSTELLATION_NAME = '32_QAM';
    Alphabet = [];
    for n = 0:1:nextpow2(2^BITS_PER_WORD)
        Alphabet = [Alphabet exp(j*pi/4).*ones(1,(nextpow2(2^5)+1))-j*n*(2*real(exp(j*pi/4))/(nextpow2(2^5)))-[(cumsum(ones(1,(nextpow2(2^5)+1)))-1).*(2*real(exp(j*pi/4))/(nextpow2(2^5)))]];
    end
    Complex_Alphabet = [Alphabet(2:5),Alphabet(7:30),Alphabet(32:35)];
    Decimal_Alphabet = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31];
elseif MODCOD == 19 %%FROM DVBT2
    BITS_PER_WORD = 4;
    CONSTELLATION_NAME = '16_QAM';
    Complex_Alphabet = [];
    for n = 0:1:(nextpow2(2^BITS_PER_WORD)-1)
        Complex_Alphabet = [Complex_Alphabet exp(j*pi/4).*ones(1,nextpow2(2^BITS_PER_WORD))-j*n*(2*real(exp(j*pi/4))/(nextpow2(2^BITS_PER_WORD)-1))-[(cumsum(ones(1,nextpow2(2^BITS_PER_WORD)))-1).*(2*real(exp(j*pi/4))/(nextpow2(2^BITS_PER_WORD)-1))]];
    end
    Complex_Alphabet = Complex_Alphabet([1 5 2 6 13 9 14 10 4 8 3 7 16 12 15 11]);
    Decimal_Alphabet = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
    Binary_Alphabet = custom_dec_to_bin(Decimal_Alphabet,BITS_PER_WORD);
end

end