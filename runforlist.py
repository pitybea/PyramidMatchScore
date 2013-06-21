import os
import sys


if __name__ == "__main__":
#    os.chdir('E:\\TunnelHough')
    
    icd=sys.argv[1]
    ifl=sys.argv[2]
    
    ins = open( ifl, "r" )
    array = []
    for line in ins:
        array.append( line.strip() )
    
    fls=array[1:int(array[0])+1]
    for fl in fls:
        cmd='c:\Python27\python.exe %s.py %s'%(icd,fl)
        print cmd
        os.system(cmd)
        
   