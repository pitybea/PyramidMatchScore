import cv2  
import sys
import os
import itertools
import numpy

def func(fl,p1,p2,p3):

        allcount_fl=open('all%s.txt'%(p2))
        allcount=int(allcount_fl.readline().strip())
        allcount_fl.close();        
        
        ptfln=open('%s_%sgtn.txt'%(fl,p1))
        pts_num=int(ptfln.readline().strip())
        ptfln.close()
        
        pts_pos=[]
        pts_pf=open('%s_%sgt.txt'%(fl,p1))
        
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
            #cv2.rectangle(img,pt1,pt2,(255,0,255),1) 
            roi=img[y if y>=0 else 0:y+40,x if x>=0 else 0:x+100]
            imn='mytest%s-%d.pgm'%(p3,allcount+i)
            cv2.imwrite(imn,roi)
            print imn
            
            
        allcount_fl=open('all%s.txt'%(p2),'w')
        allcount_fl.write(str(allcount+pts_num))
        allcount_fl.close();            
        
        #cv2.imshow("img",img)
        #cv2.waitKey(0)        

if __name__ == '__main__':
    #os.chdir('E:\\CarData\\TestImages')
    func(sys.argv[1],'','positive','pos')
    func(sys.argv[1],'neg','negative','neg')