function [mapping] = check_diploc(pos)
%--------------------------------------------------------------------------
% Till Habersetzer (25.01.2022)
% Communication Acoustics, CvO University Oldenburg
% till.habersetzer@uol.de 
%
% Dipole locations of a 2-dipol-fit are checked in terms of their mappings
% to the left and right hemisphere.
%
% Input:
%   pos: dipole locations (2,3) of 2 dipoles
%
% Output:
%   mapping: [1,2]: first dipole (1st row) is in left hemisphere
%            second dipole (2nd row) is in right hemisphere
%            [2,1]: Other way around
%   If both dipoles belong to the same hemisphere, the functions returns an
%   error.
%--------------------------------------------------------------------------
% check if source locations have different signs in x-direction

s1 = sign(pos(1,1)); % sign first position
s2 = sign(pos(2,1)); % sign second position
if ~isequal(s1,s2)
   if s1<0 %left
       mapping = [1,2];
   else
       mapping = [2,1];
   end
else
    error('Source locations are invalid. Both locations belong to same hemisphere.')
end

end