#source activate roary
#for i in $(ls novel_order/*.fna|cut -d "/" -f 2);do prokka --cdsrnaolap --cpus 10 --outdir novel_prokka/${i%%.*} --locustag ${i%%.*} --prefix ${i%%.*} novel_order/$i;done
for i in $(ls novel_prokka/);do cp novel_prokka/$i/*.gff novel_gff/;done
