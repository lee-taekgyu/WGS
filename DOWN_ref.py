 import os

l=['M']
for i in range(1,23):
    l.append(str(i))
l.append('X')
l.append('Y')
  
for i in l:
    cmd = 'wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/c    hr{}.fa.gz'.format(i)
    os.system(cmd)
for i in l:    
    cmd2 = 'cat chr{} >> hg19.fa'.format(i)
    os.system(cmd2)
