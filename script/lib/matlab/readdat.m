function [d] = readdat (fname, dim)

% [D] = READDAT(FNAME,DIM) will readin the input data in FNAME. The data
% consists of DIM-dimensional vectors.

  fid = fopen(fname, 'rt'); % Open the global input file
  if fid < 2
    fprintf(2,'[readdat] Cannot find file %s\n', fname);
    quit;
  else
  [d, counter] = fscanf(fid, '%e', [dim Inf]);
  fprintf(2, '[readdat] Read %d entries from %s\n', counter/dim, fname);
  fclose (fid);
  end

