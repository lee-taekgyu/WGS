import os

l=['MT']
for i in range(1,23):
    l.append(str(i))
l.append('X')
l.append('Y')
  
for i in l:
    cmd = 'wget ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{}.fa.gz'.format(i)
    os.system(cmd)
for i in l:    
    cmd2 = 'cat chr{} >> hg19.fa'.format(i)
    os.system(cmd2)