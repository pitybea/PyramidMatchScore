#include <stdio.h>
#include <stdlib.h>
#include <vector>
using namespace std;
#include <math.h>
#include <direct.h>
#include <Windows.h>
#include <string>
#include "../../../fileIoinclude/FileInOut.h"

int pow2(int n)
{
	return (int)pow(2.0,n);
}

struct Point
{
	double x;
	double y;
	Point(double a,double b)
	{
		x=a;
		y=b;
	}
};
struct Rect
{
	double x;
	double y;
	double width;
	double height;

	Rect(double p1,double p2,double p3,double p4)
	{
		x=p1;
		y=p2;
		width=p3;
		height=p4;
	}

	bool contains(Point pt)
	{
		return (pt.x>=x)&&(pt.x<=(x+width))&&(pt.y>=y)&&(pt.y<=(y+height));
	}
};



double overlabRatio(Rect a,Rect b)
{
	vector<Point> inps;
	inps.clear();
	vector<Point> vinps;
	vinps.clear();
	inps.push_back(Point(b.x,b.y));
	inps.push_back(Point(b.x+b.width,b.y));
	inps.push_back(Point(b.x+b.width,b.y+b.height));
	inps.push_back(Point(b.x,b.y+b.height));
	
	int situa(0);

	for (int i=0;i<inps.size();i++)
	{
		if (a.contains(inps[i]))
		{
			vinps.push_back(inps[i]);
			situa+=pow2(i);
		}
	
	}
	double ovAre(0.0);
	switch(situa)
	{
	case 0:
		if (b.contains(Point(a.x,a.y)))
		{
			ovAre=a.width*a.height;
		}
		break;
	case 1:
		ovAre=(a.x+a.width-b.x)*(a.y+a.height-b.y);
		break;
	case 2:
		ovAre=(b.x+b.width-a.x)*(a.y+a.height-b.y);
		break;
	case 3:
		ovAre=b.width*(a.y+a.height-b.y);
		break;
	case 4:
		ovAre=(b.x+b.width-a.x)*(b.y+b.height-a.y);
		break;
	
	case 6:
		ovAre=b.height*(b.x+b.width-a.x);
		break;
	case 8:
		ovAre=(a.x+a.width-b.x)*(b.y+b.height-a.y);
		break;
	case 9:
		ovAre=b.height*(a.x+a.width-b.x);
		break;
	case 12:
		ovAre=b.width*(b.y+b.height-a.y);
		break;
	case 15:
		ovAre=b.height*b.width;
		break;
	default:
		
		break;
	}
	
	return (double)ovAre/(double)(a.width*a.height+b.width*b.height-ovAre);
}

double gentheval(string s)
{

	vector<vector<double> > vals;


	vals=fileIOclass::InVectorSDouble(s+"_spvcpp.txt");


	int x,y;

	FILE* fp;
	fopen_s(&fp,(s+"_sp.txt").c_str(),"r");
	fscanf_s(fp,"%d %d",&x,&y);
	fclose(fp);

	Rect tojg((double)x,(double)y,100.0,40.0);

	double ort(-1.0);
	double finalrslt(-100);

	int thi(0);

	for (int i = 0; i < vals.size(); i++)
	{
		Rect smt( vals[i][0], vals[i][1], vals[i][2]-vals[i][0], vals[i][3]-vals[i][1]);
		double trt=overlabRatio (smt,tojg);

		if(trt>ort)
		{
			ort=trt;
			finalrslt=vals[i][5];
			thi=i;
		}

	}

	fopen_s(&fp,(s+"_vls.txt").c_str(),"w");
	fprintf(fp,"%lf %lf %lf %lf %lf %lf",ort,finalrslt,vals[thi][0],vals[thi][1],vals[thi][2],vals[thi][3]);
	fclose(fp);
	
	return finalrslt;


}

int doforlist(string s)
{
	vector<string> fnms;
	fnms.clear();
	fnms=fileIOclass::InVectorString(s.c_str());
	for(auto ss:fnms)
	{
		printf("%lf\n",	gentheval(ss));
	}

	return 0;
}

int main()
{

	_chdir("E:\\carData\\TestImages\\ps");

	doforlist("positive.lst");
	doforlist("negative.lst");

	return 0;
}