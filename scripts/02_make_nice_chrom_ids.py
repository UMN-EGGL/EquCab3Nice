
import locuspocus as lp
ids = dict()  

with open('../ensemble_id_map.txt','r') as IN:                  
    for line in IN:                                          
        i,j = line.strip().split(',')                        
        ids[i] = j     

ec3 = lp.Fasta('EquCab3')
with open('EquCab3_nice_chrom_ids.fasta','w') as OUT:                                                                                                        
    for chrom_name in ec3.chrom_names():                                        
        chrom = ec3[chrom_name]                                                 
        easy_id = ids[chrom_name]                                               
        if easy_id == 'chrUn':
            easy_id = easy_id + '_' + chrom_name
        print(f'>{easy_id} {chrom_name} {" ".join(chrom._attrs)}',file=OUT)     
        for i in range(0,len(chrom),50):                                        
            print(''.join(chrom.seq[i:i+50]),file=OUT) 

