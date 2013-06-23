#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <fstream>

#include "..\..\..\fileIoinclude\FileInOut.h"
#include <vector>
#include <string>
#include <windows.h>
#include <direct.h>
using namespace std;
using namespace cv;
//#define debugMode

class featype
{
	public:
	featype()
	{

	};
	Point pos;
	double feature[128];

	static void initOne(FILE* fp,featype &t)
	{
		fscanf(fp,"%d %d\n",&t.pos.x,&t.pos.y);
		for (int i=0;i<128;i++)
		{
			fscanf(fp,"%lf ",&t.feature[i]);
		}
	};
	static void printOne(FILE* fp, featype t)
	{
		fprintf(fp,"%d %d\n",t.pos.x,t.pos.y);
		for (int i=0;i<128;i++)
		{
			fprintf(fp,"%lf ",t.feature[i]);
		}
	};

};
int main(int argc, char* argv[])
{
	string s;
#ifdef debugMode
	_chdir("E:\\CarData\\TrainImages");
	s="pos-2";
#else
	if (argc<2)
	{
		printf("not enough parameters");
		return -1;
	}
	
	s=argv[1];
#endif	
	
	
		//string s=img_names[img_index];
		
		Mat iptimg=imread(s+".pgm");

	
		SiftFeatureDetector detector;
		vector<KeyPoint> temkpts;
		detector.detect(iptimg,temkpts);
		SiftDescriptorExtractor extractor;
		Mat descriptor;
		extractor.compute(iptimg,temkpts,descriptor);
	
	
		vector<featype> feas;
		feas.clear();
		for (auto i=0;i<descriptor.rows;i++)
		{
			featype temfes;

			temfes.pos.x=(int)temkpts[i].pt.x;
			temfes.pos.y=(int)temkpts[i].pt.y;
			Scalar ss(255,0,0);
		
			//temfes.feature.resize(128,0.0);
			double temsum(0.000001);
			for (auto j=0;j<128;j++)
			{
				temfes.feature[j]=descriptor.at<float>(i,j);
				temsum+=temfes.feature[j];
			}
			for (auto j=0;j<128;j++)
			{
				temfes.feature[j]/=temsum;
	
			}
	
	
			feas.push_back(temfes);
		}
		fileIOclass::OutVector(s+"_sift.txt",feas);
	
		
	
    return 0;

}

