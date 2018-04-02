import os 
import minus80 as m80
import locuspocus as lp

from minus80.RawFile import RawFile

ids = dict()  
with open('ensemble_id_map.txt','r') as IN:                  
    for line in IN:                                          
        i,j = line.strip().split()                        
        ids[i] = j     
if m80.Tools.available('EquCab3',dtype='Fasta'): 
    ec3 = lp.Fasta('EquCab3')
else:
    ec3 = lp.Fasta.from_file('EquCab3','EquCab3.fa.gz')

if not os.path.exists('EquCab3_nice.fasta'):
    with open('EquCab3_nice.fasta','w') as OUT:                                                                                                        
        for chrom_name in ec3.chrom_names():                                        
            print(f'Printing out {chrom_name}')
            chrom = ec3[chrom_name]                                                 
            # try to get easy_id
            if chrom_name in ids:
                chrom_name = f"{ids[chrom_name]} {chrom_name}"
            # Keep track of lengths 
            start_length = len(chrom)
            printed_length = 0
            # Loop and print
            print(f'>{chrom_name} {" ".join(chrom._attrs)}',file=OUT)     
            for i in range(0,len(chrom),70):                                        
                sequence = chrom.seq[i:i+70]
                print(''.join(sequence),file=OUT) 
                printed_length += len(sequence)
            if printed_length != start_length:
                raise ValueError('Chromosome was truncated during printing')

if not os.path.exists('EquCab3_nice.gff'):
    with RawFile('GCF_002863925.1_EquCab3.0_genomic.gff.gz') as IN, open('EquCab3_nice.gff','w') as OUT:
        for line in IN:
            id,*fields = line.split('\t')
            if id in ids:
                id = ids[id]
            print(id,*fields,file=OUT,sep='\t',end='')

