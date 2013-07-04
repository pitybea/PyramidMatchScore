import cv2  
import sys
import os
import itertools
import numpy

def func(p2):
    imagefilename = '%s.pgm'%(p2)
    img=cv2.imread(imagefilename, 1)  
    imn='%s.png'%(p2)
    cv2.imwrite(imn,img)    

if __name__ == '__main__':
    #os.chdir('E:\\CarData\\TestImages')
    func(sys.argv[1])