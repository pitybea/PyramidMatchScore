import sys
import cv2  
import sys
import os
import itertools
import numpy
from rec import *

from random import randrange

def drawrec(img,pt,sz,cl):
    y=int(pt[0])
    x=int(pt[1])
    #cv2.circle(img,p,2,(0,0,255))    
    pt1=(x,y)
    pt2=(x+sz[1],y+sz[0])
    cv2.rectangle(img,pt1,pt2,cl,1)     
if __name__ == '__main__':
    
    #for fl in fls:
    #os.chdir('E:\\carData\\TestImages')
    if(len(sys.argv)>1):
        fl=sys.argv[1]
    else:
        fl='test-99'
    
    ptfln=open('%s_gtn.txt'%(fl))
    pts_num=int(ptfln.readline().strip())
    ptfln.close()
    
    pts_pos=[]
    pts_pf=open('%s_gt.txt'%(fl))
    
    for i in range(0,pts_num):
        pts_pos.append(pts_pf.readline().strip().split(' '))
    pts_pf.close()    
    
    imagefilename = '%s.pgm'%(fl)
    img=cv2.imread(imagefilename, 1)     
    
    
    
    print img.shape
    h=img.shape[0]
    w=img.shape[1]
    
    tst=recdel()
    
    sz=(40,100)
    
    negnum=int(h*w*1.0/1.3/40/100)
    
    if(negnum<2):
        negnum=2

    psrcs=[]
    
    ngrcs=[]
    for i in range(0,pts_num):
        y=int(pts_pos[i][0])
        x=int(pts_pos[i][1])
        rc=(y,x)
        psrcs.append(rc)
        
    xrange=w-100
    yrange=h-40
    
    negcount=0
    
    while negcount<negnum:
        tx=randrange(xrange)
        ty=randrange(yrange)
        mxratio=0.0
        rc2=(ty,tx)
        for rc in psrcs:
            tratio=tst.ratio(rc,rc2,sz)
            if(tratio>mxratio):
                mxratio=tratio
        if(mxratio<0.25):
            ngrcs.append(rc2)
            negcount=negcount+1
         
        
    for rc in psrcs:
        drawrec(img,rc,sz,(255,255,0))
    for rc in ngrcs:
        drawrec(img,rc,sz,(255,0,0))
        
        
    ptfln=open('%s_neggtn.txt'%(fl),'w')
    ptfln.write(str(len(ngrcs)))
    ptfln.write('\n')    
    ptfln.close()
    
  
    pts_pf=open('%s_neggt.txt'%(fl),'w')
    for rc in ngrcs:
        pts_pf.write(' '.join(map(str,rc)))
        pts_pf.write('\n')
    pts_pf.close()
    
    #cv2.imshow("img",img)
    #cv2.waitKey(0)