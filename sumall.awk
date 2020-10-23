{
  for (i=1; i<=NF; i++) { sum[i]+= $i }   
} END { for (i=1; i<=NF; i++ ) { print "Region[", i, "] =", sum[i] } } 