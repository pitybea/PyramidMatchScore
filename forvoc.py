import cv2  
import sys
import os
import itertools
import numpy

glb=[]

def negstring(p2):
    s=(
        '# PASCAL Annotation Version 1.00\n'
        '\n'
        'Image filename : "VOC2006/PNGImages/00%d.png"\n'
        'Image size (X x Y x C) : 100 x 40 x 3\n'
        'Database : "The VOC2006 Database"\n'
        'Objects with ground truth : 0 { }\n'
        '\n'
        '# Note that there might be other objects in the image\n'
        '# for which ground truth data has not been provided.\n'
        '\n'
        '# Top left pixel co-ordinates : (1, 1)\n'%(p2)
    )
    return s

def posstring(p2):
    s=(
        '# PASCAL Annotation Version 1.00\n'
        '\n'
        'Image filename : "VOC2006/PNGImages/00%d.png"\n'
        'Image size (X x Y x C) : 100 x 40 x 3\n'
        'Database : "The VOC2006 Database"\n'
        'Objects with ground truth : 1 { "PAScar" }\n'
        '\n'
        '# Note that there might be other objects in the image\n'
        '# for which ground truth data has not been provided.\n'
        '\n'
        '# Top left pixel co-ordinates : (1, 1)\n'
        '\n'
        '# Details for object 1 ("PAScar")\n'
        'Original label for object 1 "PAScar" : "PAScar"\n'
        'Bounding box for object 1 "PAScar" (Xmin, Ymin) - (Xmax, Ymax) : (1, 1) - (99, 39)        \n'%(p2)
    )
    return s


def doone(p2,p3,inx):
    #print 'hello'
    imagefilename = '%s.pgm'%(p2)
    img=cv2.imread(imagefilename, 1) 
    cv2.imwrite('00%d.png'%(p3),img)
    fl=open('00%d.txt'%(p3),'w')
    if(inx==1):
        fl.write(posstring(p3))
    else:
        fl.write(negstring(p3))
    fl.close()
    
def dosome(p2,stn,inx):
    ifl=p2
        
    ins = open( ifl, "r" )
    array = []
    for line in ins:
        array.append( line.strip() )
    fls=array[1:int(array[0])+1]    
    
    
    i=0
    for fl in fls:
        doone(fl,stn+i,inx)
        glb.append('00%d %d\n'%(stn+i,inx))
        i=i+1
    
    

if __name__ == '__main__':
    num_stt=5304
    pos_num=550
    dosome('positive.lst',num_stt,1)
    dosome('negative.lst',num_stt+pos_num,-1)
    simfil=open('car_train.txt','w')
    for g in glb:
        simfil.write(g)
    simfil.close()