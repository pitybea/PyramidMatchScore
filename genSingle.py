import os
import sys

from rec import *

tst=recdel()
sz=(40,100)
def genpos(fl,pos,fea,opos,l1,l2,l3):
    allcount_fl=open('all%s.txt'%(l1))
    allcount=int(allcount_fl.readline().strip())
    allcount_fl.close();
    
    ptfln=open('%s_%sgtn.txt'%(fl,l3))
    pts_num=int(ptfln.readline().strip())
    ptfln.close()
    
    pts_pos=[]
    pts_pf=open('%s_%sgt.txt'%(fl,l3))
    
    for i in range(0,pts_num):
        pts_pos.append(pts_pf.readline().strip().split(' '))
    pts_pf.close()          
    
    for pti in range(0,len(pts_pos)):
        pt=pts_pos[pti]
        rc=(int(pt[0]),int(pt[1]))
        mps=[]
        mfea=[]
        for i in range(0,len(pos)):
            if(tst.ptinrec(rc,sz,pos[i])):
               mps.append('%d %d'%(pos[i][1]-rc[1],pos[i][0]-rc[0]))
               mfea.append(fea[i])
        
        mfl=open('mytest%s-%d.txt'%(l2,allcount+pti),'w')
        mfl.write(str(len(mps)))
        mfl.write('\n')
        for i in range(0,len(mps)):
            mfl.write('%s\n'%(mps[i]))
            mfl.write('%s\n'%mfea[i])
        
        mfl.close()
        mfl=open('mytest%s-%d_inx.inx'%(l2,allcount+pti),'w')  
        mfl.write('%s\n'%(fl))
        mfl.close()

    allcount_fl=open('all%s.txt'%(l1),'w')
    allcount_fl.write(str(allcount+pts_num))
    allcount_fl.close();    
    print 'ok'
    
    


def lntopt(ln):
    tm=ln.split(' ')
    return (int(tm[1]),int(tm[0]))
    

def testsub(fl):
    positions=[]
    features=[]


    temfile=open(fl+'.txt','r')
    num=int(temfile.readline().strip())
    for i in range(0, num):
        positions.append(temfile.readline().strip())
        features.append(temfile.readline().strip())
    
    pos=map(lntopt,positions)
    genpos(fl,pos,features,positions,'positive','pos','')
    genpos(fl,pos,features,positions,'negative','neg','neg')
      

if __name__ == '__main__':
    #os.chdir('E:\\CarData\\TestImages')
    testsub(sys.argv[1])
    #os.chdir('E:\\CarData\\TestImages')
    #testsub('test-26')