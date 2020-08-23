%% function that de-noises the RREF by column

function output = denoise_c(Xsm, thres)

if nargin < 2
    thres = 1e-4;
end

% Xsm: Xrref is the input RREF to be de-noised, Xsm contains only the
% non-zero rows of the RREF matrix.
% thres: the threshold parameter to be applied to each column.

output = Xsm;

% de-noising process 
for j = 1:size(Xsm, 2)
    
    output(Xsm(:,j) < thres, j) = 0;
    
end

end