
from __future__ import print_function
import os
import sys

#os.chdir('/home/pitybea/TH/originalData')
def testsub(fl):
    positions=[]
    features=[]

    temfile=open(fl+'.txt','r')
    num=int(temfile.readline().strip())
    for i in range(0, num):
        positions.append(temfile.readline().strip())
        features.append(temfile.readline().strip())
    
    temfile.close()    
    fileposs=open(fl+'_num.txt','w')
   
    print(num,file=fileposs)
    
    fileposs.close() 
    
    
    
    fileposs=open(fl+'_pos.txt','w')
   #    print(len(positions),file=fileposs)
    for pos in positions:
        print(pos,file=fileposs)
    
    fileposs.close()      
    temfile.close()    
    fileposs=open(fl+'_fea.txt','w')
   #    print(len(positions),file=fileposs)
    for fea in features:
        print(fea,file=fileposs)
    
    fileposs.close()          





if __name__ == '__main__':
    testsub(sys.argv[1])
    