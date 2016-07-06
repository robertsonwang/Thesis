% Lags variable p periods and places zeros in the rows corresponding to
% first p periods

function out = lagger(y,p)

[R,C] = size(y);
y1 = y(1:(R-p),:);
out = [zeros(p,C); y1];


