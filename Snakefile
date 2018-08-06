rule NCBI_FASTA:
	input:
		fasta=ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/latest_assembly_versions/GCF_002863925.1_EquCab3.0/GCF_002863925.1_EquCab3.0_genomic.fna.gz
	shell:
		wget {input.fasta} 

rule NCBI_GFF:
	input:
		gff=ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/latest_assembly_versions/GCF_002863925.1_EquCab3.0/GCF_002863925.1_EquCab3.0_genomic.gff.gz
	shell:
		wget {input.gff} 
