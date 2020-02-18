function nWrite = write_aeroflex_file( x, fname, verboseFlag )
%nWrite = write_aeroflex_file( x, fname, verboseFlag )
%   Scales and converts a complex, floating point file into a file consumable
%   by the Aeroflex IQCreator program to produce an AIQ file.  

if nargin < 3
  verboseFlag = false;
end

verboseFlag = verboseFlag > 0;


% Open file
if verboseFlag == true
  fprintf( '  Opening file:  ''%s'' for writing\n', fname );
end

fd = fopen( fname, 'wb' );

if fd > 0
  nWrite = fwrite( fd, format_for_aeroflex( x ), 'int16' );
  fclose( fd );
  if verboseFlag == true
    fprintf( '   Done\n' );
  end
else
  nWrite = -1;
  if verboseFlag == true
    fprintf( 2, '   Error opening file %s\n', fname );
  end
end
