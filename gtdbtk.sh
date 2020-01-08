source activate gtdbtk
gtdbtk identify --genome_dir ref_genome/ --out_dir gtdb_result/identify/ -x fna --cpus 10
gtdbtk align --identify_dir gtdb_result/identify/ --out_dir gtdb_result/align/ --cpus 10
gtdbtk infer --msa_file gtdb_result/align/gtdbtk.bac120.user_msa.fasta --out_dir gtdbtk_result/tree/ --cpus 10 
