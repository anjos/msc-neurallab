# This gawk script dries out the network file output by our training
# routines. This file contains repeatitive entries for _all_ saved networks
# during the training process. I'm only interested at the last nevertheless.

# André Rabello <Andre.Rabello@ufrj.br>

BEGIN {
  RS = "[[:space:]]*Dump";
  counter = -1;
}

{
  ++counter;
  if (counter > 0) save = $0; # The first record is not valid
}

END {
  printf "                          Dump%s", save;
}
