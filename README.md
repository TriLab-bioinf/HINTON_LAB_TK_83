## Identification of blastp-based best reciprocal matches

### Run blastp between two protein datasets
blastp -query MG_protein.fasta -db LF82_protein.fasta -evalue 1e-5 -out MG_vs_LF82.bp -outfmt '7 qaccver qlen saccver slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -mt_mode 1 -num_threads 2

### Get best reciprocal match from blastp hits with get_best_reciprocal_match.pl 
cat MG_vs_LF82.bp| ./get_best_reciprocal_match.pl > best_reciprocal_match.txt

### Adding gene names to best_reciprocal_match.txt with xlookup.pl  
xlookup.pl -d MG1655.names.txt -dri 2 -dqi 0 -i best_reciprocal_match.txt |xlookup.pl -d MG1655.names.txt -dri 3 -dqi 0| xlookup.pl -d LF82.names.txt -dri 1 -dqi 0 -ii 1 | xlookup.pl -d LF82.names.txt -dri 2 -dqi 0 -ii 1 > best_reciprocal_match_all.txt
