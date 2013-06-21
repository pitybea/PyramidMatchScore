#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <fstream>


#include <vector>
#include <string>
#include <windows.h>
#include <direct.h>
#include <string>
using namespace std;
using namespace cv;
#define debugMode


void oneImageToFeatures(string s)
{
	
	Mat iptimg=imread(s+".pgm");

	SiftFeatureDetector detector;
	vector<KeyPoint> temkpts;
	detector.detect(iptimg,temkpts);
	SiftDescriptorExtractor extractor;
	Mat descriptor;
	extractor.compute(iptimg,temkpts,descriptor);
	
	string sn,sp,sf;
	sn= s+"_num.txt";
	sp= s+"_pos.txt";
	sf= s+"_sift.txt";

	FILE* fpn=fopen(sn.c_str(),"w");
	FILE* fpp=fopen(sp.c_str(),"w");;
	FILE* fpf=fopen(sf.c_str(),"w");;
	

	fprintf(fpn,"%d\n",descriptor.rows);
	for (auto i=0;i<descriptor.rows;i++)
	{

		fprintf(fpp,"%d %d\n",(int)temkpts[i].pt.x,(int)temkpts[i].pt.y);
		//(int)temkpts[i].pt.x;
		//(int)temkpts[i].pt.y;
		Scalar ss(255,0,0);
		

		double temsum(0.000001);
		for (auto j=0;j<128;j++)
		{
			temsum+=descriptor.at<float>(i,j);
				
		}
		for (auto j=0;j<128;j++)
		{
			
			fprintf(fpf,"%lf ",descriptor.at<float>(i,j)/temsum);
		}
		fprintf(fpf,"\n");

	}
	fclose(fpn);
	fclose(fpp);
	fclose(fpf);
}

int main()
{
	_chdir("E:\\CarData\\TrainImages");	
	
	string flst;
	flst="positive.lst";
	
	FILE* fls=fopen(flst.c_str(),"r");
	int num;
	fscanf(fls,"%d\n",&num);
	
	
	printf("%d\n",num);
	for (int i = 0; i < num; i++)
	{
		//printf("%d ",i);
		char tem[50];
		fscanf(fls,"%s\n",tem);
		string ts(tem);
		printf("%s ",ts);

		try
		{
			oneImageToFeatures(ts);
		}
		catch(int m)
		{
			printf("%d ",m);
		}

		
	}
	
	fclose(fls);

	getchar();
		
	
    return 0;

}

