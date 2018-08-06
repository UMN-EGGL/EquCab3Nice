rule NCBI_FASTA:
    output:
        "data/EquCab3.ncbi.fna.gz"
    shell:
        "wget -O {output} ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/latest_assembly_versions/GCF_002863925.1_EquCab3.0/GCF_002863925.1_EquCab3.0_genomic.fna.gz"

rule NCBI_GFF:
    output:    
        "data/EquCab3.ncbi.gff.gz"
    shell:
        "wget -O {output} ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/latest_assembly_versions/GCF_002863925.1_EquCab3.0/GCF_002863925.1_EquCab3.0_genomic.gff.gz"


rule NICE_FILES:
	input:
		gff='data/EquCab3.ncbi.gff.gz',
		fna='data/EquCab3.ncbi.fna.gz'
	shell:
		"echo 1"
