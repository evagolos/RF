function [ B,ind ] = maxab( A )
% [ B,ind ] = maxab( A )
%   This function returns the largest value in data series A (or in each
%   column, if A is a matrix. This will be the largest absolute value,
%   irrespective if that is negative or positive

if all(isnan(A)), 
    B = nan; ind = [];
    return
end

if isrow(A), A = A(:); end

mabA = max(abs(A));
abA = abs(A);
B = zeros(1,size(A,2));
ind = zeros(1,size(A,2));
for ii = 1:size(A,2)
    B(ii) = unique(A(abA(:,ii)==mabA(ii),ii));
    ind(ii) = find(A(:,ii)==B(ii),1,'first');
end


end

