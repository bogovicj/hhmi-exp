function names = randomName( N, i )
% name = randomName( i )
%   i - optional number
% 
% if i is negative, it yields the number of names in the list
% 
% Grabbed the initial list from:
%   http://www.desiquintans.com/downloads/nounlist.txt

global NAMELIST;
global NUMNAMES;
    
if( ~exist('N', 'var') || isempty( N ))
    N = 1;
end

if( ~exist('i', 'var') || isempty( i ))
    if( N > 0 )
        i = randi( NUMNAMES, N, 1 );
    end
else
    N = length(i);
end

if ( N == 0 )
    [~, names] = system( sprintf( 'cat %s | wc -l', NAMELIST));
    names = str2double( names );
    return;
end

awkstr = sprintf('NR==%d||', i );
awkstr = awkstr(1:end-2);
[~, names] = system( sprintf( 'awk ''%s'' %s', awkstr, NAMELIST));
names = strsplit( strtrim( names ));