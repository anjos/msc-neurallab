function [n, x] = myhist(data,bin,area);

% [N,X] = MYHIST(DATA,BIN,AREA) will return the histogram of DATA with BIN
% bins and will normalize the area under the histogram to AREA.

  [n,x] = hist(data,bin);
  n = area*(n./sum(n));

