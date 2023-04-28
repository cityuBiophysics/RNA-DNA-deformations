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

awk  '(NR>127 && NR<147){print $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10} (NR>155 && NR<175){print $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' dnaout.lis >> curves_torsion.txt
rm dna*
done
