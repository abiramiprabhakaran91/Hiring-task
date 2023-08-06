
average_length=$(zcat NC_000913.faa.gz |awk '/^>/ {if (seqlen) { sum += seqlen; count++ } seqlen=0; next} {seqlen+=length($0)} END { sum += seqlen; count++; print sum/count }')

echo "The average length of given protein sequence = " $average_length
