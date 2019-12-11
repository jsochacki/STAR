function y = format_for_aeroflex( x )
%y = format_for_aeroflex( x )
%   Formats a complex floating-point vector into a scaled, interleaved,
%   14/16-bit signed integer vector.

nBits = 14;
scaleFactor = 2^(nBits-1)-1;

% Normalize, scale and round
mx = max(max([real(x); imag(x)]));

x = round( scaleFactor * x ./ mx );

% Convert to interleaved complex
y(2:2:2*length(x),1) = int16( imag( x ) );
y(1:2:2*length(x),1) = int16( real( x ) );
