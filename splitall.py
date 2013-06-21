import os
import sys


if __name__ == "__main__":
#    os.chdir('E:\\TunnelHough')
    
    ins = open( "positive.lst", "r" )
    array = []
    for line in ins:
        array.append( line.strip() )
    
    fls=array[1:int(array[0])+1]
    for fl in fls:
        cmd='c:\Python27\python.exe splitInfo.py %s'%(fl)
        print cmd
        os.system(cmd)
        
    ins = open( "negative.lst", "r" )
    array = []
    for line in ins:
        array.append( line.strip() )
    
    fls=array[1:int(array[0])+1]
    for fl in fls:
        cmd='c:\Python27\python.exe splitInfo.py %s'%(fl)
        print cmd
        os.system(cmd)    