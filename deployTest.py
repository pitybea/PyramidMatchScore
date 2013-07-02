import os
import sys
import glob

def printAndEx(s):
    print s
    os.system(s)
    
def initInx():
    cmd='del all*tive.txt'
    printAndEx(cmd)
    
    fl=open('allpositive.txt','wb')
    fl.write('0')
    fl.close()
    fl=open('allnegative.txt','wb')
    fl.write('0')
    fl.close()    
def deployTest():
    
    
    initInx()
    
    cmd="C:\\Python27\\python.exe runforlist.py genSingle test.lst"
    printAndEx(cmd)
    
    initInx()
        
    cmd="C:\\Python27\\python.exe runforlist.py cropImage test.lst"
    printAndEx(cmd)    
    
    cmd="del mytest\\*  /s /f  /q"
    printAndEx(cmd)    

    cmd='move mytest* mytest'
    printAndEx(cmd)
    
    cmd='copy ..\TrainImages\pca* mytest'
    printAndEx(cmd)
    
    cmd='copy *.py mytest'
    printAndEx(cmd)
    
    cmd='copy *.m mytest'
    printAndEx(cmd)
    
    os.chdir('mytest')
    
    fl=open('positive.lst','w')
    files=glob.glob("*testpos*.txt")
    fl.write('%d\n'%(len(files)))
    for fle in files:
        fl.write('%s\n'%fle.split('.')[0])
    fl.close()
    
    fl=open('negative.lst','w')
    files=glob.glob("*testneg*.txt")
    fl.write('%d\n'%(len(files)))
    for fle in files:
        fl.write('%s\n'%fle.split('.')[0])
    fl.close()    
    
    cmd='C:\\Python27\\python.exe runforlist.py splitInfo positive.lst'
    printAndEx(cmd)
    
    cmd='C:\\Python27\\python.exe runforlist.py splitInfo negative.lst'
    printAndEx(cmd)
    

if __name__ == '__main__':
    deployTest()