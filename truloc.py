import os
import sys

if __name__ == '__main__':
    os.chdir('E:\\CarData')
    temfile=open('trueLocations.txt','r')
    
    lns=[]
    for i in range(0,170):
        temln=temfile.readline().strip()
       
        lns.append([x.strip() for x in temln.replace(' ','').replace(')(',':').replace('(','').replace(')','').split(':')])
    
    temfile.close()
    for ln in lns:
        fl=open('test-'+ln[0]+'_gt.txt','w')
        fl1=open('test-'+ln[0]+'_gtn.txt','w')
        fl1.write(str(len(ln)-1))
        if(len(ln)>1):
            for j in range(1,len(ln)):
                fl.write(ln[j].replace(',',' '))
                fl.write('\n')
        
        fl.close()
        fl1.close()
        
    #print lns
   