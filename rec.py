import os
import sys

class recdel:
    
    def ptinrec(self,rc,sz,pt):
        x=pt[1]
        y=pt[0]
        
        h=sz[0]
        w=sz[1]
        
        rx=rc[1]
        ry=rc[0]
        return (x>=rx) & (x<=(rx+w)) & (y>=ry) & (y<=(ry+h))
    
    def ratio(self,rc1,rc2,sz):
        y1=rc1[0]
        x1=rc1[1]
        y2=rc2[0]
        x2=rc2[1]
        
        h=sz[0]
        w=sz[1]
        
        pts=[]
        pts.append((y1,x1))
        pts.append((y1,x1+w))
        
        pts.append((y1+h,x1+w))
        pts.append((y1+h,x1))
        
        are=2*w*h
        
        tem =-1

        #print pts
        
        for i in range(0,len(pts)):
            if(self.ptinrec(rc2,sz,pts[i])):
                #print i
                tem=i
                break
                
        if (tem==-1):
            return 0.0
        elif (tem==0):
            ov = (y2+h-y1)*(x2+w-x1)*1.0
        elif (tem==1):
            ov= (y2+h-y1)*(x1+w-x2)*1.0
        elif(tem==2):
            ov=(y1+h-y2)*(x1+w-x2)*1.0
        elif(tem==3):
            ov= (y2+h-y1)*(x2+w-x1)*1.0
        
        #print ov,are
        return ov*1.0/(are-ov)


if __name__ == '__main__':
    a=recdel()
    print a.ratio((10,10),(10,10),(20,10))    
    
    