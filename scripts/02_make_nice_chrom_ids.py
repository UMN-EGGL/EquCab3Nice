
import minus80 as m80
import locuspocus as lp

ids = dict()  
with open('ensemble_id_map.txt','r') as IN:                  
    for line in IN:                                          
        i,j = line.strip().split(',')                        
        ids[i] = j     
if m80.Tools.available('EquCab3',dtype='Fasta'): 
    ec3 = lp.Fasta('EquCab3')
else:
    ec3 = lp.Fasta.from_file('EquCab3','EquCab3.fa.gz')

with open('EquCab3_nice_chrom_ids.fasta','w') as OUT:                                                                                                        
    for chrom_name in ec3.chrom_names():                                        
        chrom = ec3[chrom_name]                                                 
        easy_id = ids[chrom_name]                                               
        if easy_id == 'chrUn':
            easy_id = easy_id + '_' + chrom_name
        print(f'>{easy_id} {chrom_name} {" ".join(chrom._attrs)}',file=OUT)     
        for i in range(0,len(chrom),70):                                        
            print(''.join(chrom.seq[i:i+70]),file=OUT) 

