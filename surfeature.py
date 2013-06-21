import sys


import cv2  
import sys
import os
import numpy
#os.chdir('E:\\TunnelHough')

def main(argv):
    if len(argv) < 2:
        fl='testdata\\frame92'
    else:
        fl=argv[1]
    img=cv2.imread(fl+'.jpg', 1)  
    hessian_threshold = 60
    surf = cv2.SURF(hessian_threshold)
    imgg =cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)
    kp, descritors = surf.detect(imgg,None,useProvidedKeypoints = False)
    
    
        #cv2.circle(img,(int(x),int(y)),2,(255,255,0))
        
    rowsize = len(descritors) / len(kp)
    if rowsize > 1:
        rows = numpy.array(descritors, dtype = numpy.float32).reshape((-1, rowsize))      
    print rows.shape
    fil=open('%s_num.txt'%fl,'wb')
    fil.write('%d'%(len(kp)))
    fil.close()
    
    fil=open('%s_pos.txt'%fl,'wb')
    for i in range(0,len(kp)):
        x,y=kp[i].pt
        fil.write('%d %d\n'%(x,y))
    fil.close()    

    fil=open('%s_fea.txt'%fl,'wb')
    for row in rows:
        fil.write(' '.join(map(str, row)) )
        fil.write('\n')
    fil.close()    
    
if __name__ == "__main__":
    main(sys.argv)
