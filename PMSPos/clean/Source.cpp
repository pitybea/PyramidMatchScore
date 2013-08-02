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


void newgivescoreshelp(PMStruc pedmd,int dim,string s)
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
		vector<double> result;
		double score;

		vector<vector<double> > vecs;
		vecs=selectVecButLstTwo ( tvd,dim);
		vector<double> reslt;
		//pedmd.givePyramidMatchScore(vecs,false,reslt);
		
		auto st=toTwo(vecs,2);

		vector<vector<double> > app=st.second;
		vector<vector<double> >ps=st.first;
	

		vector<vector<pair<int,int> > > fss;
		for (int i = 0; i < app.size(); i++)
		{
			fss.push_back(pedmd.dataToTwoInxs(app[i]));
		}

		score=pedmd.givePyramidMatchScore(fss,ps,false,reslt);
	//	if(!another)
		//	score=pedmd.givePyramidMatchScore(selectVecButLstTwo ( tvd,dim),false,result);
	
		
		printf("%lf\n",score);
	
	}
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
		vector<double> result;
		double score;
	
		score=pedmd.givePyramidMatchScore(selectVecButLstTwo ( tvd,dim),false,result);
	
	
		printf("%lf\n",score);

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
void givescoressplit(PMStruc pedmd,int dim,bool ada,bool another)
{
	
	
	
	adaboost machine;

//	FILE* fp=fopen("adabMach.txt","r");
//	machine.loadFromfile(fp);
//	fclose(fp);
	newgivescoreshelp(pedmd,dim,"positive.lst");
	newgivescoreshelp(pedmd,dim,"negative.lst");
	

	
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
		PMStruc pedmd(PMStruc::postitionSpecific);
		printf("-------------******************-----(%d)--------**********************\n",dim);
		pedmd.generatePymFromdata(selectVec(data,dim),dim-2);
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


	return 0;
}
int test_example_num()
{
	_chdir("E:\\carData\\TrainImages");
	
	auto ad=getAllfeas(10);
	auto data=ad.first;


	PMSEnsemble pem;
	pem.generateAaBsFromdata(data,6);

	PMStruc n(PMStruc::postitionSpecific);
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

	pmse.generateAaBsFromdata(data,6);
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

	

	PMStruc pedmd(PMStruc::postitionSpecific);
	printf("-------------******************-----(%d)--------**********************\n",dim);
	pedmd.generatePymFromdata(data,dim);
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



	PMStruc pedmd(PMStruc::postitionSpecific);

	pedmd.generatePymFromdata(data,dim);
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



	PMStruc pedmd(PMStruc::postitionSpecific);

	pedmd.generatePymFromdata(selectVec(data,dim),dim);
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
	pedmd.generatePymFromdata(data,dim);
	vector<string> flNms;
	flNms.clear();
	flNms=fileIOclass::InVectorString("positive.lst");
	
	for(auto s : flNms)
	{
		vector<vector<double> > tvd;
		tvd.clear();
		tvd=fileIOclass::InVectorSDouble(s+"_ptscpp.txt");
		double siz=(double)tvd.size();
		vector<double> result;
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
		vector<double> result;
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
	
	PMStruc pedmd(PMStruc::postitionSpecific);
	//printf("-------------******************-----(%d)--------**********************\n",dim);
	pedmd.generatePymFromdata(data,dim);
	vector<string> flNms;
	flNms.clear();
	flNms=fileIOclass::InVectorString("positive.lst");
	
	for(auto s : flNms)
	{
		vector<vector<double> > tvd;
		tvd.clear();
		tvd=fileIOclass::InVectorSDouble(s+"_ptscpp.txt");
		double siz=(double)tvd.size();
		vector<double> result;
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
		vector<double> result;
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
	pedmd.generatePymFromdata(data,dim);
	givescores(pedmd,dim,false,false);

	return 0;
}

int benchmark(int dim)
{
	_chdir("E:\\carData\\TrainImages");
	
//	int dim=10;
	auto ad=getAllfeas(dim);
	auto data=ad.first;
	
	_chdir("E:\\carData\\TestImages\\mytest");



	PMStruc pedmd(PMStruc::positionsimple,5);
	printf("-------------******************-----(%d)--------**********************\n",dim);
	pedmd.generatePymFromdata(data,dim);
	givescores(pedmd,dim,false,false);

	return 0;
}

int splitmethod()
{
	_chdir("E:\\carData\\TrainImages");
	
	int dim=10;
	auto ad=getAllfeas(dim);
	auto data=ad.first;
	
	_chdir("E:\\carData\\TestImages\\mytest");



	PMStruc pedmd(PMStruc::positionsimple,5);
//	printf("-------------******************-----(%d)--------**********************\n",dim);
	pedmd.generatePymFromdata(data,dim);
	givescoressplit(pedmd,dim,false,false);

	return 0;
}

int testposLvlLmt()
{
	_chdir("E:\\carData\\TrainImages");
	
		int dim=10;
	auto ad=getAllfeas(dim);
	auto data=ad.first;
	
	_chdir("E:\\carData\\TestImages\\mytest");


	for (int i = 3; i < 8; i++)
	{
		PMStruc pedmd(PMStruc::postitionSpecific,i);
		
		pedmd.generatePymFromdata(data,dim);
		printf("-------------******************-----(%d)--------**********************\n",i);
		givescores(pedmd,dim,false,false);
	}
	

	return 0;
}





int testposspe(PMStruc pedmd,int dim,double ratio)
{
	pedmd.GeneratePosWeightWithParameter(ratio);
	givescores(pedmd,dim,false,false);

	return 0;
}

int testPosWeight()
{

	_chdir("E:\\carData\\TrainImages");
	
	int dim=10;
	auto ad=getAllfeas(dim);
	auto data=ad.first;
	
	_chdir("E:\\carData\\TestImages\\mytest");



	PMStruc pedmd(PMStruc::postitionSpecific);
	pedmd.generatePymFromdata(data,dim);


	double ra=0.5;
	
	for (int i = 0; i < 20; i++)
	{
		ra=0.05*(i);
		printf("-------------******************-----(%lf)--------**********************\n",ra);

		testposspe(pedmd,dim,ra);
	}
	
	return 0;
}

int testposLvlLmtRatio()
{
	_chdir("E:\\carData\\TrainImages");
	
		int dim=10;
	auto ad=getAllfeas(dim);
	auto data=ad.first;
	
	_chdir("E:\\carData\\TestImages\\mytest");

	double ra=0.5;
	for (int i = 3; i < 8; i++)
	{
		PMStruc pedmd(PMStruc::postitionSpecific,i);
		
		pedmd.generatePymFromdata(data,dim);
		printf("-------------******************-----(%d)--------**********************\n",i);
		for (int j = 0; j < 20; j++)
		{
			ra=0.05*(j);
			printf("-------------******************-----(%lf)--------**********************\n",ra);

			testposspe(pedmd,dim,ra);
		}
	
	}
	

	return 0;
}

int formatchange(string s,int pnum)
{
	vector<vector<double> > allvals;
	allvals=fileIOclass::InVectorSDouble(s.c_str());

	for (int i = 0; i < allvals.size(); i++)
	{
		if (i<pnum)
		{
			printf("+1 ");
		}
		else
			printf("-1 ");
		for (int j = 0; j < allvals[i].size(); j++)
		{
			printf("%d:%lf ",j+1,allvals[i][j]);
		}
		printf("\n");
	}
	return 0;
}

int trainweights()
{
	_chdir("E:\\project\\PyramidMatchScore\\PMS\\x64\\Release");
	vector<vector<double> > allvals;
	allvals=fileIOclass::InVectorSDouble("training.txt");
	vector<vector<double> > pvals;
	pvals.clear();
	vector<vector<double> >nvals;
	nvals.clear();

	pvals.insert(pvals.end(),allvals.begin(),allvals.begin()+550);
	nvals.insert(nvals.end(),allvals.begin()+550,allvals.end());

	vector<vector<double> > lhbds;
	vector<double> stps;
	lhbds.resize(2,vector<double> (allvals[0].size(),0.0));
	stps.resize(allvals[0].size(),0.0);

	for (int i = 0; i < allvals[0].size(); i++)
	{
		double tm=pow(2.0,i+1);
		lhbds[0][i]=0-tm;
		lhbds[1][i]=tm;
		stps[i]=tm/30.0;
	}

	trainVecs vec(trainVecs::allAtOnce,lhbds,stps,pvals,nvals);
	vec.train();
	printf("-------------******************-----()--------**********************\n");
	vec.printWeights();
	return 0;
}

int main()
{

//	vector<vector<pair<int,int> > > r;
//	r=genPAorders(6);

//	PMStruc pedmd(PMStruc::postitionSpecific,5);
//	benchmark();


	benchmark(1);
	benchmark(3);
	benchmark(7);
	benchmark(10);
	benchmark(15);
	benchmark(20);
//	splitmethod();
//	test_dimension_bynumber();
	return 0;
}