import sys
import cv2  
import sys
import os
import itertools
import numpy

if __name__ == '__main__':
    
    #for fl in fls:
    fl=sys.argv[1]
    
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
    for i in range(0,pts_num):
        y=int(pts_pos[i][0])
        x=int(pts_pos[i][1])
        #cv2.circle(img,p,2,(0,0,255))    
        pt1=(x,y)
        pt2=(x+100,y+40)
        cv2.rectangle(img,pt1,pt2,(255,0,255),1) 
        
        
    ptfln=open('%s_neggtn.txt'%(fl))
    pts_num=int(ptfln.readline().strip())
    ptfln.close()
    
    pts_pos=[]
    pts_pf=open('%s_neggt.txt'%(fl))
    
    for i in range(0,pts_num):
        pts_pos.append(pts_pf.readline().strip().split(' '))
    pts_pf.close()    
    
   
    for i in range(0,pts_num):
        y=int(pts_pos[i][0])
        x=int(pts_pos[i][1])
        #cv2.circle(img,p,2,(0,0,255))    
        pt1=(x,y)
        pt2=(x+100,y+40)
        cv2.rectangle(img,pt1,pt2,(255,0,0),1)         
     
    
    
        
    ptfln=open('%s_num.txt'%(fl))
    pts_num=int(ptfln.readline().strip())
    ptfln.close()
    
    pts_pos=[]
    pts_pf=open('%s_pos.txt'%(fl))
    
    for i in range(0,pts_num):
        pts_pos.append(pts_pf.readline().strip().split(' '))
    pts_pf.close()    
    
   
    for i in range(0,pts_num):
        p=(int(pts_pos[i][0]),int(pts_pos[i][1]))
        cv2.circle(img,p,2,(0,0,255))    
        
    cv2.imshow("img",img)
    cv2.waitKey(0)