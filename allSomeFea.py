

from __future__ import print_function
import os
import sys
if __name__ == "__main__":
#    os.chdir('E:\\TunnelHough')
    
    ins = open( "negative.lst", "r" )
    array = []
    for line in ins:
        array.append( line.strip() )
    
    fls=array[1:int(array[0])+1]
    features=[]
    positions=[]
    for fl in fls:             
        temfile=open(fl+'.txt','r')
        num=int(temfile.readline().strip())
        for i in range(0, num):
            positions.append(temfile.readline().strip())
            features.append(temfile.readline().strip())
        
        temfile.close()   
    
    allf=open('allNegativeFeatures.txt','wb')
    for fea in features:
        print(fea,file=allf)
    
    allf.close()
        