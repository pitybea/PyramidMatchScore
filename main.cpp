/*
 * main.cpp
 *
 *  Created on: Jul 13, 2012
 *      Author: pitybea
 */


#include "../../../fileIoinclude/FileInOut.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <math.h>
#include <stdio.h>
#include <Windows.h>
#include <direct.h>
#include "PMS.h"







void getMatric(vector<vector<int> >&vls, string fname)
{
	ifstream iof;
	iof.open(fname.c_str());
	int wid,hei;
	iof>>wid>>hei;
	vls.clear();
	vls.resize(hei,vector<int>(wid,0));
	for (int i=0;i<hei;i++)
	{
		for (int j=0;j<wid;j++)
		{
			iof>>vls[i][j];
		}

	}

}
void ptsTov(vector<Point> p,vector<vector<double> >& d)
{
	d.resize(p.size(),vector<double>(2,0.0));
	for (int i=0;i<p.size();i++)
	{
		d[i][0]=p[i].x;
		d[i][1]=p[i].y;
	}

}

int main()
{
	_chdir("E:\\uiucCars\\CarData\\TrainImages");

	vector<featype> allfeas;
	allfeas.clear();
	allfeas=fileIOclass::InVector<featype>("allPostiveFeatures.txt");
//	fileIOclass::CombineFromFileList<featype>("neg.txt","sift.txt","allNegtiveFeatures.txt");

	return 0;

}


