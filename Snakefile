import gzip
import bz2
import lzma
import os 

samples = [
"11F_ATGTCA_L007_R1_001","11F_ATGTCA_L007_R2_001", "11F_ATGTCA_L008_R1_001", "11F_ATGTCA_L008_R2_001",
"11M_GGCTAC_L004_R1_001","11M_GGCTAC_L004_R2_001", "12F_GTCCGC_L007_R1_001", "12F_GTCCGC_L007_R2_001",
"12F_GTCCGC_L008_R1_001","12F_GTCCGC_L008_R2_001", "12M_AGTCAA_L005_R1_001", "12M_AGTCAA_L005_R2_001",
"12M_AGTCAA_L006_R1_001","12M_AGTCAA_L006_R2_001", "14F_GATCAG_L005_R1_001", "14F_GATCAG_L005_R2_001",
"14F_GATCAG_L006_R1_001","14F_GATCAG_L006_R2_001", "14M_ACAGTG_L004_R1_001", "14M_ACAGTG_L004_R2_001",
"17F_CCGTCC_L007_R1_001","17F_CCGTCC_L007_R2_001", "17F_CCGTCC_L008_R1_001", "17F_CCGTCC_L008_R2_001",
"17M_CTTGTA_L005_R1_001","17M_CTTGTA_L005_R2_001", "17M_CTTGTA_L006_R1_001", "17M_CTTGTA_L006_R2_001",
"1F_TAGCTT_L007_R1_001","1F_TAGCTT_L007_R2_001", "1F_TAGCTT_L008_R1_001", "1F_TAGCTT_L008_R2_001",
"1M_GCCAAT_L004_R1_001","1M_GCCAAT_L004_R2_001", "20F_CAGATC_L005_R1_001", "20F_CAGATC_L005_R2_001",
"20F_CAGATC_L006_R1_001","20F_CAGATC_L006_R2_001", "20M_TTAGGC_L004_R1_001", "20M_TTAGGC_L004_R2_001",
"35F_CTTGTA_L007_R1_001","35F_CTTGTA_L007_R2_001", "35F_CTTGTA_L008_R1_001", "35F_CTTGTA_L008_R2_001",
"35M_ACTTGA_L004_R1_001","35M_ACTTGA_L004_R2_001", "36F_GTGGCC_L007_R1_001", "36F_GTGGCC_L007_R2_001",
"36F_GTGGCC_L008_R1_001","36F_GTGGCC_L008_R2_001", "36M_ATGTCA_L005_R1_001", "36M_ATGTCA_L005_R2_001",
"36M_ATGTCA_L006_R1_001","36M_ATGTCA_L006_R2_001", "37F_AGTCAA_L007_R1_001", "37F_AGTCAA_L007_R2_001",
"37F_AGTCAA_L008_R1_001","37F_AGTCAA_L008_R2_001", "37M_GATCAG_L004_R1_001", "37M_GATCAG_L004_R2_001",
"38F_TGACCA_L007_R1_001","38F_TGACCA_L007_R2_001", "38F_TGACCA_L008_R1_001", "38F_TGACCA_L008_R2_001",
"38M_ATTCCT_L005_R1_001","38M_ATTCCT_L005_R2_001", "38M_ATTCCT_L006_R1_001", "38M_ATTCCT_L006_R2_001",
"39F_ATCACG_L007_R1_001","39F_ATCACG_L007_R2_001", "39F_ATCACG_L008_R1_001", "39F_ATCACG_L008_R2_001",
"39M_CGTACG_L005_R1_001","39M_CGTACG_L005_R2_001", "39M_CGTACG_L006_R1_001", "39M_CGTACG_L006_R2_001",
"40F_CGTACG_L007_R1_001","40F_CGTACG_L007_R2_001", "40F_CGTACG_L008_R1_001", "40F_CGTACG_L008_R2_001",
"40M_GTCCGC_L005_R1_001","40M_GTCCGC_L005_R2_001", "40M_GTCCGC_L006_R1_001", "40M_GTCCGC_L006_R2_001",
"41F_ACTTGA_L005_R1_001","41F_ACTTGA_L005_R2_001", "41F_ACTTGA_L006_R1_001", "41F_ACTTGA_L006_R2_001",
"41M_TGACCA_L004_R1_001","41M_TGACCA_L004_R2_001", "49F_ACAGTG_L007_R1_001", "49F_ACAGTG_L007_R2_001",
"49F_ACAGTG_L008_R1_001","49F_ACAGTG_L008_R2_001", "49M_ATCACG_L005_R1_001", "49M_ATCACG_L005_R2_001",
"49M_ATCACG_L006_R1_001","49M_ATCACG_L006_R2_001", "50F_GCCAAT_L005_R1_001", "50F_GCCAAT_L005_R2_001",
"50F_GCCAAT_L006_R1_001","50F_GCCAAT_L006_R2_001", "50M_CGATGT_L004_R1_001", "50M_CGATGT_L004_R2_001",
"51F_GAGTGG_L007_R1_001","51F_GAGTGG_L007_R2_001", "51F_GAGTGG_L008_R1_001", "51F_GAGTGG_L008_R2_001",
"51M_GTGAAA_L005_R1_001","51M_GTGAAA_L005_R2_001", "51M_GTGAAA_L006_R1_001", "51M_GTGAAA_L006_R2_001",
"52F_CGATGT_L007_R1_001","52F_CGATGT_L007_R2_001", "52F_CGATGT_L008_R1_001", "52F_CGATGT_L008_R2_001",
"52M_GAGTGG_L005_R1_001","52M_GAGTGG_L005_R2_001", "52M_GAGTGG_L006_R1_001", "52M_GAGTGG_L006_R2_001",
"61F_TTAGGC_L007_R1_001","61F_TTAGGC_L007_R2_001", "61F_TTAGGC_L008_R1_001", "61F_TTAGGC_L008_R2_001",
"61M_ACTGAT_L005_R1_001","61M_ACTGAT_L005_R2_001", "61M_ACTGAT_L006_R1_001", "61M_ACTGAT_L006_R2_001",
"64F_ACTGAT_L007_R1_001","64F_ACTGAT_L007_R2_001", "64F_ACTGAT_L008_R1_001", "64F_ACTGAT_L008_R2_001",
"64M_GTGGCC_L005_R1_001","64M_GTGGCC_L005_R2_001", "64M_GTGGCC_L006_R1_001", "64M_GTGGCC_L006_R2_001",
"65F_ATTCCT_L007_R1_001","65F_ATTCCT_L007_R2_001", "65F_ATTCCT_L008_R1_001", "65F_ATTCCT_L008_R2_001",
"65M_GTTTCG_L005_R1_001","65M_GTTTCG_L005_R2_001", "65M_GTTTCG_L006_R1_001", "65M_GTTTCG_L006_R2_001",
"66F_AGTTCC_L007_R1_001","66F_AGTTCC_L007_R2_001", "66F_AGTTCC_L008_R1_001", "66F_AGTTCC_L008_R2_001",
"66M_TAGCTT_L004_R1_001","66M_TAGCTT_L004_R2_001", "67F_GTGAAA_L007_R1_001", "67F_GTGAAA_L007_R2_001",
"67F_GTGAAA_L008_R1_001","67F_GTGAAA_L008_R2_001", "67M_AGTTCC_L005_R1_001", "67M_AGTTCC_L005_R2_001",
"67M_AGTTCC_L006_R1_001","67M_AGTTCC_L006_R2_001", "68F_GTTTCG_L007_R1_001", "68F_GTTTCG_L007_R2_001",
"68F_GTTTCG_L008_R1_001","68F_GTTTCG_L008_R2_001", "68M_CCGTCC_L005_R1_001", "68M_CCGTCC_L005_R2_001",
"68M_CCGTCC_L006_R1_001","68M_CCGTCC_L006_R2_001", "69F_CAGATC_L007_R1_001", "69F_CAGATC_L007_R2_001",
"69F_CAGATC_L008_R1_001","69F_CAGATC_L008_R2_001", "69M_TTAGGC_L005_R1_001", "69M_TTAGGC_L005_R2_001",
"69M_TTAGGC_L006_R1_001","69M_TTAGGC_L006_R2_001", "82F_ACTTGA_L007_R1_001", "82F_ACTTGA_L007_R2_001",
"82F_ACTTGA_L008_R1_001","82F_ACTTGA_L008_R2_001", "82M_TGACCA_L005_R1_001", "82M_TGACCA_L005_R2_001",
"82M_TGACCA_L006_R1_001","82M_TGACCA_L006_R2_001", "86F_GCCAAT_L007_R1_001", "86F_GCCAAT_L007_R2_001",
"86F_GCCAAT_L008_R1_001","86F_GCCAAT_L008_R2_001", "86M_CGATGT_L005_R1_001", "86M_CGATGT_L005_R2_001",
"86M_CGATGT_L006_R1_001","86M_CGATGT_L006_R2_001", "87F_GGCTAC_L007_R1_001", "87F_GGCTAC_L007_R2_001",
"87F_GGCTAC_L008_R1_001","87F_GGCTAC_L008_R2_001", "87M_CAGATC_L004_R1_001", "87M_CAGATC_L004_R2_001",
"90F_ACAGTG_L005_R1_001","90F_ACAGTG_L005_R2_001", "90F_ACAGTG_L006_R1_001", "90F_ACAGTG_L006_R2_001",
"90M_ATCACG_L004_R1_001","90M_ATCACG_L004_R2_001" 
]

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

rule all:
    input:
        S3.remote(expand('HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}.fastq.gz',sample=samples))
        

rule compress_fastq:
    input:
        S3.remote("RnaSeqData/Project_McCue_Project_022/{sample}.fastq")
    output:
        S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}.fastq.gz')
    run:
        import ipdb; ipdb.set_trace()
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

