import gzip
import bz2
import lzma
import os 



'''
The file hierarchy on S3 looks like:

    s3://HorseGeneAnnotation/
        private/
            sequence/
                RNASEQ/
                    fastq/
                    bam/
                WGS/
                    fastq/
                    bam/
            variant/
                vcf/
        public/
            refgen/
                GCF_002863925.1_EquCab3.0/
                    GCF_002863925.1_EquCab3.0.nice.fna
                    GCF_002863925.1_EquCab3.0.nice.nice.gff
                    GCF_002863925.1_EquCab3.0.nice.index
                        
'''

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
    'NC_009144.3':'chr1',  'NC_009145.3':'chr2',  'NC_009146.3':'chr3',
    'NC_009147.3':'chr4',  'NC_009148.3':'chr5',  'NC_009149.3':'chr6',
    'NC_009150.3':'chr7',  'NC_009151.3':'chr8',  'NC_009152.3':'chr9',
    'NC_009153.3':'chr10', 'NC_009154.3':'chr11', 'NC_009155.3':'chr12',
    'NC_009156.3':'chr13', 'NC_009157.3':'chr14', 'NC_009158.3':'chr15',
    'NC_009159.3':'chr16', 'NC_009160.3':'chr17', 'NC_009161.3':'chr18',
    'NC_009162.3':'chr19', 'NC_009163.3':'chr20', 'NC_009164.3':'chr21',
    'NC_009165.3':'chr22', 'NC_009166.3':'chr23', 'NC_009167.3':'chr24',
    'NC_009168.3':'chr25', 'NC_009169.3':'chr26', 'NC_009170.3':'chr27',
    'NC_009171.3':'chr28', 'NC_009172.3':'chr29', 'NC_009173.3':'chr30',
    'NC_009174.3':'chr31', 'NC_009175.3':'chrX',  'NC_001640.1':'chrMt'
}


GCF = 'GCF_002863925.1_EquCab3.0'


projects,samples, = S3.glob_wildcards("RnaSeqData/{project}/{sample}.fastq")
sample_dict = {s:p for p,s in zip(projects,samples)}


rule all:
    input:
        S3.remote(expand('HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}.fastq.gz',sample=samples))

rule compress_fastq:
    input:
        [S3.remote(f"RnaSeqData/{x}/{y}.fastq") for x,y in zip(projects,samples)]
    output:
        S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}.fastq.gz')
    run:
        shell('gzip -c {input} > {output}')

rule STAR_INDEX:
    input:
        gff = S3.remote("HorseGeneAnnotation/public/refgen/{GCF}/{GCF}_genomic.nice.gff.gz"),
        fna = S3.remote("HorseGeneAnnotation/public/refgen/{GCF}/{GCF}_genomic.nice.fna.gz")
    output:
        S3.remote('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/Genome'),
        S3.remote('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/SA'),
        S3.remote('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/SAindex'),
        S3.remote('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/chrLength.txt'),
        S3.remote('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/chrName.txt'),
        S3.remote('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/chrNameLength.txt'),
        S3.remote('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/chrStart.txt'),
        S3.remote('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/exonGeTrInfo.tab'),
        S3.remote('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/exonInfo.tab'),
        S3.remote('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/geneInfo.tab'),
        S3.remote('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/genomeParameters.txt'),
        S3.remote('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/sjdbInfo.txt'),
        S3.remote('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/sjdbList.fromGTF.out.tab'),
        S3.remote('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/sjdbList.out.tab'),
        S3.remote('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/transcriptInfo.tab') 
    threads: 30
    shell:
        '''STAR \
          --runThreadN {threads} \
          --runMode genomeGenerate \
          --genomeDir HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/ \
          --genomeFastaFiles {input.fna} \
          --sjdbGTFfile {input.gff} \
          --sjdbGTFtagExonParentTranscript Parent 
        '''

rule NICE_FASTA:
    input:
        fna = FTP.remote('ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/all_assembly_versions/{GCF}/{GCF}_genomic.fna.gz')
    output:
        fna = S3.remote("HorseGeneAnnotation/public/refgen/{GCF}/{GCF}_genomic.nice.fna.gz")
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
        gff = FTP.remote('ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/all_assembly_versions/{GCF}/{GCF}_genomic.gff.gz')
    output:
        gff = S3.remote("HorseGeneAnnotation/public/refgen/{GCF}/{GCF}_genomic.nice.gff.gz")
    run:
        with RawFile(input.gff) as IN, \
            open(output.gff,'w') as OUT:
            for line in IN:
                id,*fields = line.split('\t')
                if id in id_map:
                    id = id_map[id]
                print(id,*fields,file=OUT,sep='\t',end='')

