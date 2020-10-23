#!/bin/bash

#echo -n "Enter the starting line number: "
#read s
#echo -n "Enter the end line number: "
#read e
echo -n "Enter the number of lines from the end: "
read number

tail -n $number lightcurve.dat > temp.dat
s=1
e=$number
gawk -v start=$s -v fin=$e '
{ 
if(NR>=start && NR<=fin)
{ 
  for (i = 1; i <= NF; i++) 
  {
    array[i]+=$i;
    if(NR == start) {
    startime=$1;
                    }
  if(NR == fin) {
    endtime=$1;
                }
   }
}
}
END { 
for(x=2; x<=NF; x++)
  {
  print array[x]/(endtime-startime);
  }
print "start time [s] = ", startime;
print "end time [s] = ", endtime;
print "span [s] =", endtime-startime;
}
' temp.dat > t.dat
#sed 'N;$!P;$!D;$d' t.dat > t1.dat
##delete the last three lines
sed -n -e :a -e '1,3!{P;N;D;};N;ba' t.dat > t1.dat
paste nuRange.dat t1.dat > nuVsFA.dat
tail -n 3 t.dat
rm temp.dat
rm t.dat
rm t1.dat
