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
#include "adaboost.h"
#include "../../../fileIoinclude/FileInOut.h"




pair<vector<vector<double> >,vector<vector<vector<double> > > > getAllfeas(int dim)
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

		auto tvm=selectVecButLstTwo(tvd,dim);

		rslt.insert(rslt.end(),tvm.begin(),tvm.end());
		arst.push_back(tvm);
	}


	return pair<vector<vector<double> >,vector<vector<vector<double> > > >(rslt,arst);
	
}

void givescoreshelp(PMStruc pedmd,int dim,string s,bool ada,bool another,adaboost machine)
{
	vector<string> flNms;
	flNms.clear();
	flNms=fileIOclass::InVectorString(s.c_str());
	for(auto s : flNms)
	{
		vector<vector<double> > tvd;
		
		tvd.clear();
		tvd=fileIOclass::InVectorSDouble(s+"_ptscpp.txt");
		double siz=(double)tvd.size();
		vector<int> result;
		double score;
		if(!another)
			score=pedmd.givePyramidMatchScore(selectVecButLstTwo ( tvd,dim),false,result);
		else
			score=pedmd.givePyramidMatchScore(selectVec ( tvd,dim),false,result);
		if(ada)
		{
			vector<double> to;
			to.resize(result.size(),0.0);
			for (int i = 0; i < result.size(); i++)
			{
				to[i]=(double)result[i]/siz;
				printf("%lf ",to[i]);
			}
			printf("%lf\n",machine.classfy(to));
		}
		else
		{
			vector<double> to;
			to.resize(result.size(),0.0);
			for (int i = 0; i < result.size(); i++)
			{
				to[i]=(double)result[i]/siz;
				printf("%lf ",to[i]);
			}
			printf("%lf\n",score);
		}
	}
}

void givescores(PMStruc pedmd,int dim,bool ada,bool another)
{
	
	
	
	adaboost machine;

//	FILE* fp=fopen("adabMach.txt","r");
//	machine.loadFromfile(fp);
//	fclose(fp);
	givescoreshelp(pedmd,dim,"positive.lst",ada,another,machine);
	givescoreshelp(pedmd,dim,"negative.lst",ada,another,machine);
	

	
}



void givescores(PMSEnsemble pedmd,int dim)
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
		printf("%lf\n",pedmd.givePyramidMatchScore(selectVecButLstTwo( tvd,dim)));
		
	}

	flNms.clear();
	flNms=fileIOclass::InVectorString("negative.lst");
	
	for(auto s : flNms)
	{
		vector<vector<double> > tvd;
		tvd.clear();
		tvd=fileIOclass::InVectorSDouble(s+"_ptscpp.txt");
		vector<int> result;
		printf("%lf\n",pedmd.givePyramidMatchScore(selectVecButLstTwo( tvd,dim)));
		
	}
}


void testDimension(pair<vector<vector<double> >,vector<vector<vector<double> > > > ad,int dim)
{
		auto data=ad.first;
		PMStruc pedmd(PMStruc::normal);
		printf("-------------******************-----(%d)--------**********************\n",dim);
		pedmd.generatePymFromdata(selectVecButLstTwo(data,dim));
		givescores(pedmd,dim,false,false);
}

int main1()
{
	vector<vector<int> > tst;
	vector<int> ok;
	ok.resize(10,0);
	for (int i = 0; i < 10; i++)
	{
		ok[i]=i;
	}

	for (int i = 0; i < 10; i++)
	{
		tst.push_back(ok);
	}

	for(auto ss:tst)
	{
		for(auto s:ss)
			printf("%d ",s);
		printf("\n");
	}
	for(auto ss:selectVecButLstTwo(tst,3))
	{
		for(auto s:ss)
			printf("%d ",s);
		printf("\n");
	}

	getchar();

	return 0;
}

int test_dimension_bynumber()
{
	_chdir("E:\\carData\\TrainImages");
	
	auto ad=getAllfeas(32);
	auto data=ad.first;
	
	_chdir("E:\\carData\\TestImages\\mytest");
	testDimension(ad,1);

	testDimension(ad,3);

	testDimension(ad,7);
	testDimension(ad,10);
	testDimension(ad,15);
	testDimension(ad,20);
	testDimension(ad,30);

	return 0;
}
int test_example_num()
{
	_chdir("E:\\carData\\TrainImages");
	
	auto ad=getAllfeas(10);
	auto data=ad.first;


	PMSEnsemble pem;
	pem.generateAaBsFromdata(data);

	PMStruc n(PMStruc::normal);
	n.initPymWithABs(pem.aAbs,data[0].size());


	_chdir("E:\\carData\\TestImages\\mytest");

	int dim=10;

	printf("-------------******************-------------**********************\n");
	for (int i = 0; i < 100; i++)
	{
		n.AddSeverlData(selectVecButLstTwo(ad.second[i],dim),true);
	}
	givescores(n,dim,false,false);

	for (int i =100; i < 200; i++)
	{
		n.AddSeverlData(selectVecButLstTwo(ad.second[i],dim),true);
	}
	printf("-------------******************-------------**********************\n");
	givescores(n,dim,false,false);



	for (int i =200; i < 300; i++)
	{
		n.AddSeverlData(selectVecButLstTwo(ad.second[i],dim),true);
	}
	printf("-------------******************-------------**********************\n");
	givescores(n,dim,false,false);


	
	for (int i =300; i < 400; i++)
	{
		n.AddSeverlData(selectVecButLstTwo(ad.second[i],dim),true);
	}
	printf("-------------******************-------------**********************\n");
	givescores(n,dim,false,false);


	
	for (int i =400; i < ad.second.size(); i++)
	{
		n.AddSeverlData(selectVecButLstTwo(ad.second[i],dim),true);
	}
	printf("-------------******************-------------**********************\n");
	givescores(n,dim,false,false);



	return 0;
}
int test_ensemble_bynumber()
{

	_chdir("E:\\carData\\TrainImages");
	
	int dim=10;
	auto ad=getAllfeas(dim);
	auto data=ad.first;

	PMSEnsemble pmse;
	pmse.threshold=2.5;

	pmse.generateAaBsFromdata(data);
	pmse.generateStructureFromData(ad.second);

	

	_chdir("E:\\carData\\TestImages\\mytest");
	givescores(pmse,dim);
	printf("-------------******************-------------**********************\n");
	printf("%d\n",pmse.pyms.size());
	for (auto pm: pmse.pyms)
	{
		printf("%d\n",pm.getNumofData());
	}

	printf("-------------******************-------------**********************\n");

		pmse.threshold=0.4;

	pmse.pyms.clear();
	pmse.generateStructureFromData(ad.second);
	givescores(pmse,dim);
	printf("-------------******************-------------**********************\n");
	printf("%d\n",pmse.pyms.size());
	for (auto pm: pmse.pyms)
	{
		printf("%d\n",pm.getNumofData());
	}


	printf("-------------******************-------------**********************\n");

			pmse.threshold=0.3;

	pmse.pyms.clear();
	pmse.generateStructureFromData(ad.second);
	givescores(pmse,dim);
	printf("-------------******************-------------**********************\n");
	printf("%d\n",pmse.pyms.size());
	for (auto pm: pmse.pyms)
	{
		printf("%d\n",pm.getNumofData());
	}

	

	//pedmd.outToAFile("pospym.txt");


	//PMStruc ptem;
	//ptem.loadFromAFile("pospym.txt");

	return 0;
}



int test_basic()
{
	_chdir("E:\\carData\\TrainImages");
	
	int dim=10;
	auto ad=getAllfeas(dim);
	auto data=ad.first;
	
	_chdir("E:\\carData\\TestImages\\mytest");

	

	PMStruc pedmd(PMStruc::normal);
	printf("-------------******************-----(%d)--------**********************\n",dim);
	pedmd.generatePymFromdata(data);
	givescores(pedmd,dim,false,false);

	return 0;
}

int no_app()
{
	
	_chdir("E:\\carData\\TrainImages");
	
		int dim=0;
	auto ad=getAllfeas(dim);
	auto data=ad.first;
	
	_chdir("E:\\carData\\TestImages\\mytest");



	PMStruc pedmd(PMStruc::normal);

	pedmd.generatePymFromdata(data);
	givescores(pedmd,dim,false,false);
	
	return 0;
}

int no_pos()
{
	
	_chdir("E:\\carData\\TrainImages");
	
		int dim=10;
	auto ad=getAllfeas(dim);
	auto data=ad.first;
	
	_chdir("E:\\carData\\TestImages\\mytest");



	PMStruc pedmd(PMStruc::normal);

	pedmd.generatePymFromdata(selectVec(data,dim));
	givescores(pedmd,dim,false,true);
	
	return 0;
}

int genpostraining()
{
	_chdir("E:\\carData\\TrainImages");
	
	int dim=10;
	auto ad=getAllfeas(dim);
	auto data=ad.first;
	
	PMStruc pedmd(PMStruc::postitionSpecific);
	//printf("-------------******************-----(%d)--------**********************\n",dim);
	pedmd.generatePymFromdata(data);
	vector<string> flNms;
	flNms.clear();
	flNms=fileIOclass::InVectorString("positive.lst");
	
	for(auto s : flNms)
	{
		vector<vector<double> > tvd;
		tvd.clear();
		tvd=fileIOclass::InVectorSDouble(s+"_ptscpp.txt");
		double siz=(double)tvd.size();
		vector<int> result;
		pedmd.givePyramidMatchScore(selectVecButLstTwo ( tvd,dim),true,result);
		for(auto ss:result)
			printf("%lf ",(double)ss/siz);
		printf("\n");
		
	}

	flNms.clear();
	flNms=fileIOclass::InVectorString("negative.lst");
	
	for(auto s : flNms)
	{
		vector<vector<double> > tvd;
		tvd.clear();
		tvd=fileIOclass::InVectorSDouble(s+"_ptscpp.txt");
		double siz=(double)tvd.size();
		vector<int> result;
		pedmd.givePyramidMatchScore(selectVecButLstTwo( tvd,dim),false,result);
		for(auto ss:result)
			printf("%lf ",(double)ss/siz);
		printf("\n");
		
	}



	return 0;
}


int gentraining()
{
	_chdir("E:\\carData\\TrainImages");
	
	int dim=10;
	auto ad=getAllfeas(dim);
	auto data=ad.first;
	
	PMStruc pedmd(PMStruc::normal);
	//printf("-------------******************-----(%d)--------**********************\n",dim);
	pedmd.generatePymFromdata(data);
	vector<string> flNms;
	flNms.clear();
	flNms=fileIOclass::InVectorString("positive.lst");
	
	for(auto s : flNms)
	{
		vector<vector<double> > tvd;
		tvd.clear();
		tvd=fileIOclass::InVectorSDouble(s+"_ptscpp.txt");
		double siz=(double)tvd.size();
		vector<int> result;
		pedmd.givePyramidMatchScore(selectVecButLstTwo( tvd,dim),true,result);
		for(auto ss:result)
			printf("%lf ",(double)ss/siz);
		printf("\n");
		
	}

	flNms.clear();
	flNms=fileIOclass::InVectorString("negative.lst");
	
	for(auto s : flNms)
	{
		vector<vector<double> > tvd;
		tvd.clear();
		tvd=fileIOclass::InVectorSDouble(s+"_ptscpp.txt");
		double siz=(double)tvd.size();
		vector<int> result;
		pedmd.givePyramidMatchScore(selectVecButLstTwo( tvd,dim),false,result);
		for(auto ss:result)
			printf("%lf ",(double)ss/siz);
		printf("\n");
		
	}



	return 0;
}

int testposspe()
{
	_chdir("E:\\carData\\TrainImages");
	
		int dim=10;
	auto ad=getAllfeas(dim);
	auto data=ad.first;
	
	_chdir("E:\\carData\\TestImages\\mytest");



	PMStruc pedmd(PMStruc::postitionSpecific);
//	printf("-------------******************-----(%d)--------**********************\n",dim);
	pedmd.generatePymFromdata(data);
	givescores(pedmd,dim,false,false);

	return 0;
}
int main()
{
	no_app();
	return 0;
}