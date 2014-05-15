function [file, dirdata] = sortdirbynum(s, f, sep)
% SORTDIRBYNUM  List the files in a directory ordering them by a numerical
% substring field in the file name
%
% [FILE, DIRDATA] = SORTDIRBYNUM(S, F, SEP)
%
%    FILE is a listing with the file names sorted accordingly.
%
%    DIRDATA is a string with the path to the target directory.
%
%    S is a string with the file names to list. Pathnames and wildcards may
%    be used, the same as with function DIR().
%
%    F is the index of the field used to sort the file names. The first
%    field has index 1. By default, F=1. The field will be treated as a
%    number, not a string. Thus, '1.' will precede '10.'.
%
%    SEP is the character that separates fields. By default, SEP='.'.
%
%    This function will also detect whether the field numeration is not
%    consecutive, or whether there are duplicates, and issue a warning.
%
%    Example:
%
%    >> SORTDIRBYNUM('MISAS_DT_22_10_2009.MR.HEAD_GENERAL.6.*.IMA', 5, '.')
%
%    will sort files 
%
%    MISAS_DT_22_10_2009.MR.HEAD_GENERAL.6.10.2009.10.23.16.02.38.562500.9699412.BMP
%    MISAS_DT_22_10_2009.MR.HEAD_GENERAL.6.11.2009.10.23.16.02.38.562500.9699451.BMP
%    MISAS_DT_22_10_2009.MR.HEAD_GENERAL.6.1.2009.10.23.16.02.38.562500.9699056.BMP
%    MISAS_DT_22_10_2009.MR.HEAD_GENERAL.6.12.2009.10.23.16.02.38.562500.9699490.BMP
%
%    as
%
%    MISAS_DT_22_10_2009.MR.HEAD_GENERAL.6.1.2009.10.23.16.02.38.562500.9699056.BMP
%    MISAS_DT_22_10_2009.MR.HEAD_GENERAL.6.10.2009.10.23.16.02.38.562500.9699412.BMP
%    MISAS_DT_22_10_2009.MR.HEAD_GENERAL.6.11.2009.10.23.16.02.38.562500.9699451.BMP
%    MISAS_DT_22_10_2009.MR.HEAD_GENERAL.6.12.2009.10.23.16.02.38.562500.9699490.BMP

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2009 University of Oxford
% Version: 0.1.0
% $Rev$
% $Date$
% 
% University of Oxford means the Chancellor, Masters and Scholars of
% the University of Oxford, having an administrative office at
% Wellington Square, Oxford OX1 2JD, UK. 
%
% This file is part of Gerardus.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. The offer of this
% program under the terms of the License is subject to the License
% being interpreted in accordance with English Law and subject to any
% action against the University of Oxford being under the jurisdiction
% of the English Courts.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


% check arguments
error( nargchk( 1, 3, nargin, 'struct' ) );
error( nargoutchk( 0, 2, nargout, 'struct' ) );

% defaults
if ( nargin < 2 || isempty( f ) )
    f = 1;
end
if ( nargin < 3 || isempty( sep ) )
    sep = '.';
end

% get path to the data files
[ dirdata, name ] = fileparts( s );

% get list of image files
file = dir( s );

if isempty(file)
    return
end

% create vector for the field that we are going to use to sort the file
% names
num = zeros( length( file ), 1 );

% extract the numerical value of the field under consideration
for I = 1:length( file )
    
    % get starting indices of any occurence of the separator
    idx = strfind( file( I ).name, sep );
    
    % add starting index so that the first substring is the first field
    idx = [ 1 idx ];
    
    % extract the desired field
    subs = file( I ).name( idx(f)+1:idx(f+1)-1 );
    
    % convert to numeric (error if non-numeric)
    aux = str2double( subs );
    if isempty( aux )
        error( [ 'Selected field ' num2str( f ) ...
            ' is non-numeric in file ' file( I ).name ] )
    end
    num( I ) = aux;
    
end

% sort the selected field in ascending order (we cannot assume that it's
% going to be from 1 to N)
[ num_sorted, idx ] = sort( num, 1, 'ascend' );

% check whether there are any slices missing
if ( length( file ) < num_sorted(end) - num_sorted(1) + 1 )
    warning( 'GERARDUS:fieldgap', 'Note that the values in the selected field have 1 or more gaps' )
end

% check whether there are any repeated field values
if ( any( diff( num_sorted ) == 0 ) )
    warning( 'GERARDUS:fielddup', 'Note that there are files with repeated values in the selected field' )
end

% sort filenames in increasing order according to the selected field
file = file( idx );
