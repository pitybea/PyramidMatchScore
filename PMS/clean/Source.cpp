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
#include "../../../fileIoinclude/FileInOut.h"

pair<vector<vector<double> >,vector<vector<vector<double> > > > getAllfeas()
{
	vector<vector<double> > rslt;
	rslt.clear();

	vector<vector<vector<double> > > arst;
	arst.clear();

	vector<string> flNms;
	flNms.clear();
	flNms=fileIOclass::InVectorString("positive.lst");
	
	for(auto s : flNms)
	{
		vector<vector<double> > tvd;
		tvd.clear();
		tvd=fileIOclass::InVectorSDouble(s+"_ptscpp.txt");

		rslt.insert(rslt.end(),tvd.begin(),tvd.end());
		arst.push_back(tvd);
	}


	return pair<vector<vector<double> >,vector<vector<vector<double> > > >(rslt,arst);
	
}




void givescores(PMStruc pedmd)
{
	
	vector<string> flNms;
	flNms.clear();
	flNms=fileIOclass::InVectorString("positive.lst");
	
	for(auto s : flNms)
	{
		vector<vector<double> > tvd;
		tvd.clear();
		tvd=fileIOclass::InVectorSDouble(s+"_ptscpp.txt");
		vector<int> result;
		printf("%lf\n",pedmd.givePyramidMatchScore(tvd,false,result));
		
	}

	flNms.clear();
	flNms=fileIOclass::InVectorString("negative.lst");
	
	for(auto s : flNms)
	{
		vector<vector<double> > tvd;
		tvd.clear();
		tvd=fileIOclass::InVectorSDouble(s+"_ptscpp.txt");
		vector<int> result;
		printf("%lf\n",pedmd.givePyramidMatchScore(tvd,false,result));
		
	}
}


void givescores(PMSEnsemble pedmd)
{
	
	vector<string> flNms;
	flNms.clear();
	flNms=fileIOclass::InVectorString("positive.lst");
	
	for(auto s : flNms)
	{
		vector<vector<double> > tvd;
		tvd.clear();
		tvd=fileIOclass::InVectorSDouble(s+"_ptscpp.txt");
		vector<int> result;
		printf("%lf\n",pedmd.givePyramidMatchScore(tvd));
		
	}

	flNms.clear();
	flNms=fileIOclass::InVectorString("negative.lst");
	
	for(auto s : flNms)
	{
		vector<vector<double> > tvd;
		tvd.clear();
		tvd=fileIOclass::InVectorSDouble(s+"_ptscpp.txt");
		vector<int> result;
		printf("%lf\n",pedmd.givePyramidMatchScore(tvd));
		
	}
}

int main ()
{
	_chdir("E:\\carData\\TrainImages");
	
	auto ad=getAllfeas();
	auto data=ad.first;


	PMSEnsemble pem;
	pem.generateAaBsFromdata(data);

	PMStruc n(PMStruc::normal);
	n.initPymWithABs(pem.aAbs,data[0].size());


	_chdir("E:\\carData\\TestImages\\mytest");

	printf("-------------******************-------------**********************\n");
	for (int i = 0; i < 100; i++)
	{
		n.AddSeverlData(ad.second[i],true);
	}
	givescores(n);

	for (int i =100; i < 200; i++)
	{
		n.AddSeverlData(ad.second[i],true);
	}
	printf("-------------******************-------------**********************\n");
	givescores(n);



	for (int i =200; i < 300; i++)
	{
		n.AddSeverlData(ad.second[i],true);
	}
	printf("-------------******************-------------**********************\n");
	givescores(n);


	
	for (int i =300; i < 400; i++)
	{
		n.AddSeverlData(ad.second[i],true);
	}
	printf("-------------******************-------------**********************\n");
	givescores(n);


	
	for (int i =400; i < ad.second.size(); i++)
	{
		n.AddSeverlData(ad.second[i],true);
	}
	printf("-------------******************-------------**********************\n");
	givescores(n);



	return 0;
}
int main_ok()
{

	_chdir("E:\\carData\\TrainImages");
	
	auto ad=getAllfeas();
	auto data=ad.first;

	PMSEnsemble pmse;
	pmse.threshold=15;

	pmse.generateAaBsFromdata(data);
	pmse.generateStructureFromData(ad.second);

	


	/*
	PMSEnsemble pem;
	pem.generateAaBsFromdata(data);



	PMStruc pedmd(PMStruc::normal);
	pedmd.generatePymFromdata(data);


		

	_chdir("E:\\carData\\TestImages\\mytest");
	givescores(pedmd);

	printf("-------------******************-------------**********************\n");


	
	*/

	_chdir("E:\\carData\\TestImages\\mytest");
	givescores(pmse);
	printf("-------------******************-------------**********************\n");
	printf("%d\n",pmse.pyms.size());
	for (auto pm: pmse.pyms)
	{
		printf("%d\n",pm.getNumofData());
	}

	printf("-------------******************-------------**********************\n");

		pmse.threshold=40;

	pmse.pyms.clear();
	pmse.generateStructureFromData(ad.second);
	givescores(pmse);
	printf("-------------******************-------------**********************\n");
	printf("%d\n",pmse.pyms.size());
	for (auto pm: pmse.pyms)
	{
		printf("%d\n",pm.getNumofData());
	}


	printf("-------------******************-------------**********************\n");

			pmse.threshold=100;

	pmse.pyms.clear();
	pmse.generateStructureFromData(ad.second);
	givescores(pmse);
	printf("-------------******************-------------**********************\n");
	printf("%d\n",pmse.pyms.size());
	for (auto pm: pmse.pyms)
	{
		printf("%d\n",pm.getNumofData());
	}

	getchar();

	//pedmd.outToAFile("pospym.txt");


	//PMStruc ptem;
	//ptem.loadFromAFile("pospym.txt");

	return 0;
}