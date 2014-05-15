function [sol] = ls_solaat(rhs)

% Yin Zhang, 1997

global probData

% ----- check input rhs -----
if issparse(rhs)
   error('RHS vector must be in full matrix format.');
end;
m = length(probData.PERM);
if (max(size(rhs)) ~= m) || (min(size(rhs)) ~= 1)
   error('No symbolic factor or input sizes mismatch.');
end;

% ------ back solve ------
sol = zeros(m,1);
sol(probData.PERM) = ls_blkslv(probData.XLNZ,probData.XSUPER,probData.XLINDX,probData.LINDX,probData.LNZ0,rhs(probData.PERM));
