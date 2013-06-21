#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <windows.h>
#include <direct.h>
using namespace std;
using namespace cv;

string Which_type(string a,vector<string> tps)
{
	for each(string s in tps)
	{
		if (a.find(s)<a.length())
		{
			return s;
		}
	}
	return string("");
}
struct featype
{
	Point pos;
	string type;
	string origImg;
	vector<float> feature;
};
int main(int argc, char* argv[])
{
	
	_chdir("D:\\share\\CarData\\my\\TrainImages");
	ifstream inpF;
	inpF.open("tI.txt");
	int n_imgs;
	inpF>>n_imgs;
	vector<string> img_names;
	img_names.resize(n_imgs,"");
	vector<string> img_labels;
	img_labels.resize(n_imgs,"");
	for (auto i=0;i<n_imgs;i++)
	{
		inpF>>img_names[i]>>img_labels[i];
	}
	inpF.close();
	
	vector<featype> codes;
	codes.clear();
	
	ofstream outF;
	outF.open("Siftforppp.txt");
//	outF<<img_names.size()<<endl;

	for (int img_index=0;img_index<img_names.size();img_index++)
	{
		string s=img_names[img_index];
		
		Mat iptimg=imread(s);

		int imgw=iptimg.cols;
		imgw/=2;
		
		int imgh=iptimg.rows;
		imgh/=2;

		SiftFeatureDetector detector;
		vector<KeyPoint> temkpts;
		detector.detect(iptimg,temkpts);
		SiftDescriptorExtractor extractor;
		Mat descriptor;
		extractor.compute(iptimg,temkpts,descriptor);
		cout<<s<<" "<<descriptor.rows<<endl;
	//	outF<<descriptor.rows<<" "<<img_labels[img_index]<<" "<<s<<" "<<imgw<<" "<<imgh<<endl;
		for (auto i=0;i<descriptor.rows;i++)
		{
			featype temfes;
			temfes.type=img_labels[img_index];
			temfes.pos.x=imgw-(int)temkpts[i].pt.x;
			temfes.pos.y=imgh-(int)temkpts[i].pt.y;

		//	outF<<temfes.pos.x<<" "<<temfes.pos.y<<endl;
			Scalar ss(255,0,0);
		
			temfes.feature.resize(128,0.0);
			float temsum(0.0);
			for (auto j=0;j<128;j++)
			{
				temfes.feature[j]=descriptor.at<float>(i,j);
				temsum+=temfes.feature[j];
			}
			for (auto j=0;j<128;j++)
			{
				temfes.feature[j]/=temsum;
				outF<<temfes.feature[j]<<" ";
			}
			outF<<endl;
			temfes.origImg=s;
			
		}
		
	}
	
	
	
	outF.close();
    return 0;

}

