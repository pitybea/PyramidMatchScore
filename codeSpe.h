#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;
using namespace cv;

struct feainfo
{

	int id;
	int imgid;
	Point pos;
	vector<double> feacontent;
};
struct imginfo
{
	int imgid;
	string imgname;
	string lable;
	int imgwidth;
	int imgheight;
};
int getFeainfo(vector<feainfo>& finfo,vector<imginfo>& iminfo)
{

	;
	finfo.clear();
	;
	iminfo.clear();

	ifstream feaInfof,imgInfof;
	feaInfof.open("feature_info.txt");
	imgInfof.open("img_info.txt");
	//pFeaf.open("pureFeature.txt");
	//feaInfof<<finfo.size()<<endl;
	int temi;
	feaInfof>>temi;
	feainfo temfi;
	finfo.resize(temi,temfi);

	for (int i=0;i<finfo.size();i++)
	{
		//feainfo temfinfo=finfo[i];
		int id;
		feaInfof>>id;
		feaInfof>>finfo[id].imgid>>finfo[id].pos.x>>finfo[id].pos.y;
		finfo[id].feacontent.resize(128,0.0);
		for (int j=0;j<finfo[id].feacontent.size();j++)
		{
			feaInfof>>finfo[id].feacontent[j];


		}
	}

	imgInfof>>temi;
	imginfo temimi;

	iminfo.resize(temi,temimi);
	//		<<iminfo.size()<<endl;
	for (int i=0;i<iminfo.size();i++)
	{

		//	temimfinfo=iminfo[i];
		int temid;
		imgInfof>>temid;
		imgInfof>>iminfo[temid].imgname>>iminfo[temid].lable>>iminfo[temid].imgwidth>>iminfo[temid].imgheight;
	}
	feaInfof.close();
	imgInfof.close();

	return 0;
}

