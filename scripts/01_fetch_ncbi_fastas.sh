parallel wget ftp://ftp.ncbi.nih.gov/genomes/Equus_caballus/CHR_0{}/eca_ref_EquCab3.0_chr{}.fa.gz ::: $(seq 1 10)
parallel wget ftp://ftp.ncbi.nih.gov/genomes/Equus_caballus/CHR_{}/eca_ref_EquCab3.0_chr{}.fa.gz ::: $(seq 11 31)
parallel wget ftp://ftp.ncbi.nih.gov/genomes/Equus_caballus/CHR_{}/eca_ref_EquCab3.0_chr{}.fa.gz ::: MT UN X
parallel wget ftp://ftp.ncbi.nih.gov/genomes/Equus_caballus/CHR_{}/eca_ref_EquCab3.0_chr{}.fa.gz ::: 10
parallel wget ftp://ftp.ncbi.nih.gov/genomes/Equus_caballus/CHR_{}/eca_ref_EquCab3.0_chr{}.fa.gz ::: Un
zcat $(ls -1 *.fa.gz | sort -V)  | gzip  > EquCab3.fa
