% Gauss elimination for tridiagonlaized system
function [cnew] = tridiag(diag1, diag2, diag3, rs, npx);

% Reduce to two-diagnoal system
for i= 2:npx,
   diag2(i) = diag2(i) - diag3(i-1)*diag1(i)/diag2(i-1);
   rs(i) = rs(i) - rs(i-1)*diag1(i)/diag2(i-1);
end

% Backward substitution
rs(npx) = rs(npx)/diag2(npx);
for i=npx-1: -1: 1,
   rs(i) = (rs(i) - diag3(i)*rs(i+1))/diag2(i);
end   

cnew = rs;

