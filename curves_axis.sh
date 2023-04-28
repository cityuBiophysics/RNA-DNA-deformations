for i in `seq 1 60000`; do
  echo $i
  cp curves_twist_pdb/dnaMD$i.pdb dna.pdb
./../../../../../curves+/Cur+ <<!
&inp file=dna, lis=dnaout,
 lib=../../../../../curves+/standard, &end
2 1 -1 0 0
1:25
50:26
!

awk  '(NR>34 && NR<53){print $5 "\t" $6 "\t" $7 "\t" $8}' dnaout.lis >> curves_axis.txt
rm dna*
done
