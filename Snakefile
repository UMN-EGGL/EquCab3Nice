import gzip
import bz2
import lzma
import os 

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
s3_key_id = os.environ.get('AWS_ACCESS_KEY')
s3_access_key = os.environ.get('AWS_SECRET_KEY')

S3 = S3RemoteProvider(
    endpoint_url='https://s3.msi.umn.edu',
    access_key_id=s3_key_id, 
    secret_access_key=s3_access_key
)

id_map = {
    'NC_009144.3':'chr1',
    'NC_009145.3':'chr2',
    'NC_009146.3':'chr3',
    'NC_009147.3':'chr4',
    'NC_009148.3':'chr5',
    'NC_009149.3':'chr6',
    'NC_009150.3':'chr7',
    'NC_009151.3':'chr8',
    'NC_009152.3':'chr9',
    'NC_009153.3':'chr10',
    'NC_009154.3':'chr11',
    'NC_009155.3':'chr12',
    'NC_009156.3':'chr13',
    'NC_009157.3':'chr14',
    'NC_009158.3':'chr15',
    'NC_009159.3':'chr16',
    'NC_009160.3':'chr17',
    'NC_009161.3':'chr18',
    'NC_009162.3':'chr19',
    'NC_009163.3':'chr20',
    'NC_009164.3':'chr21',
    'NC_009165.3':'chr22',
    'NC_009166.3':'chr23',
    'NC_009167.3':'chr24',
    'NC_009168.3':'chr25',
    'NC_009169.3':'chr26',
    'NC_009170.3':'chr27',
    'NC_009171.3':'chr28',
    'NC_009172.3':'chr29',
    'NC_009173.3':'chr30',
    'NC_009174.3':'chr31',
    'NC_009175.3':'chrX',
    'NC_001640.1':'chrMt'
}

class RawFile(object):
    def __init__(self,filename):
        self.filename = filename
        if filename.endswith('.gz'):
            self.handle = gzip.open(filename,'rt')
        elif filename.endswith('bz2'):
            self.handle = bz2.open(filename,'rt')
        elif filename.endswith('xz'):
            self.handle = lzma.open(filenaem,'rt')
        else:
            self.handle = open(filename,'r')
    def __enter__(self):
        return self.handle
    def __exit__(self,dtype,value,traceback):
        self.handle.close()


rule NICE_FASTA:
    input:
        fna='data/EquCab3.ncbi.fna.gz',
    output:
        fna='data/EquCab3.nice.fna',
    run:
        with RawFile(input.fna) as IN, open(output.fna,'w') as OUT:                                                                                                        
            for line in IN:                                        
                if line.startswith('>'):
                    name, *fields = line.lstrip('>').split()
                    if name in id_map:
                        new_name = '>' + id_map[name]
                        line = ' '.join([new_name, name] + fields + ['\n'])
                print(line,file=OUT,end='')

rule NICE_GFF:
    input:
        gff='data/EquCab3.ncbi.gff.gz'
    output:
        gff='data/EquCab3.nice.gff'
    run:

        with RawFile(input.gff) as IN, \
            open(output.gff,'w') as OUT:
            for line in IN:
                id,*fields = line.split('\t')
                if id in id_map:
                    id = id_map[id]
                print(id,*fields,file=OUT,sep='\t',end='')

rule NCBI_FASTA:
    input:
        FTP.remote("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/latest_assembly_versions/GCF_002863925.1_EquCab3.0/GCF_002863925.1_EquCab3.0_genomic.fna.gz")
    output:
        "data/EquCab3.ncbi.fna.gz"
    shell:
        'cp {input} {output}'

rule NCBI_GFF:
    input:
        FTP.remote("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/latest_assembly_versions/GCF_002863925.1_EquCab3.0/GCF_002863925.1_EquCab3.0_genomic.gff.gz")
    output:    
        "data/EquCab3.ncbi.gff.gz"
    shell:
        'cp {input} {output}'

rule UPLOAD_GFF:
    input:
        gff='data/EquCab3.nice.gff'
    output:
        S3.remote('GFFs/EquCab3.nice.gff')        
    shell:
        'cp {input} {output}'

rule UPLOAD_FASTA:
    input:
        fna='data/EquCab3.nice.fna'
    output:
        S3.remote('FASTAs/EquCab3.nice.fna')
    shell:
        'cp {input} {output}'


rule STAR_INDEX:
    input:
        fasta='data/EquCab3.nice.fna',
        gff='data/EquCab3.nice.gff'
    output:

    shell:
        '''STAR \
          --runThreadN 2 \
          --runMode genomeGenerate \
          --genomeDir /home/schae234/Codes/BuildNiceEquCab3Fasta/data/ \
          --genomeFastaFiles /home/schae234/Codes/BuildNiceEquCab3Fasta/{input.fasta} \
          --sjdbGTFfile /home/schae234/Codes/BuildNiceEquCab3Fasta/{input.gff} \
          --sjdbGTFtagExonParentTranscript Parent \
        '''

