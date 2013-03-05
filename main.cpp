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


static vector<int> oneNumber(string s,PMStruc ptem,vector<vector<double> > TrmTx)
{
	string ts=s+"sift.txt";
	vector<featype> vf=fileIOclass::InVector<featype> (ts);
	vector<vector<double> > tomatch=prepareData(vf,TrmTx);
	vector<int> result;
	ptem.givePyramidMatchScore(tomatch,true,result);

	return result;
}

void trainPositive()
{
	vector<featype> allfeas;
	allfeas.clear();
	allfeas=fileIOclass::InVector<featype>("allPostiveFeatures.txt");
	
	vector<vector<double> > TrmTx;
	TrmTx=fileIOclass::InVectorSDouble("../DimensionReduction.txt");

	TrmTx=keepFirSevDims(TrmTx,10);
	vector<vector<double > >  dats;
	dats.resize(allfeas.size(),vector<double>(0,0.0));
	for (int i=0;i<dats.size();i++)
	{
		dats[i]=allfeas[i].toVdouble();
	}
	
	ZeroMnVec(dats);
	dats=TransitMtx(dats,TrmTx);
	
	dats=addPositionsToData(dats,allfeas);
	//NormalVec(dats);
	

	PMStruc pedmd(5,PMStruc::normal,"pospym");
	pedmd.generatePymFromdata(dats);
	pedmd.outToAFile("pospym.txt");
}
void trainNegative()
{
	vector<featype> allfeas;
	allfeas.clear();
	allfeas=fileIOclass::InVector<featype>("allNegtiveFeatures.txt");

	vector<vector<double> > TrmTx;
	TrmTx=fileIOclass::InVectorSDouble("../DimensionReduction.txt");

	TrmTx=keepFirSevDims(TrmTx,10);
	vector<vector<double > >  dats;
	dats.resize(allfeas.size(),vector<double>(0,0.0));
	for (int i=0;i<dats.size();i++)
	{
		dats[i]=allfeas[i].toVdouble();
	}
	
	ZeroMnVec(dats);
	dats=TransitMtx(dats,TrmTx);
	
	//dats=addPositionsToData(dats,allfeas);
	//NormalVec(dats);
	//fileIOclass::OutVectorSDouble("temp.txt",dats,true);

	PMStruc pedmd(5,PMStruc::normal,"negpym");
	pedmd.generatePymFromdata(dats);
	pedmd.outToAFile("negpym.txt");
}




int main()
{
	_chdir("E:\\uiucCars\\CarData\\TrainImages");
//	trainPositive();
//	trainNegative();

	vector<vector<double> > TrmTx;
	TrmTx=fileIOclass::InVectorSDouble("../DimensionReduction.txt");
	TrmTx=keepFirSevDims(TrmTx,10);


	PMStruc ptem;
	ptem.loadFromAFile("pospym.txt");

	vector<string> posInames=fileIOclass::InVectorString("pos.txt");
	vector<string> negInames=fileIOclass::InVectorString("neg.txt");
	
	vector<vector<int> > posnumbers;
	vector<vector<int> > negnumbers;
	for(string s:posInames)
		posnumbers.push_back(oneNumber(s,ptem,TrmTx));

	fileIOclass::OutVectorSInt("temposnums",posnumbers,true);

	for(string s:negInames)
		negnumbers.push_back(oneNumber(s,ptem,TrmTx));
	
	fileIOclass::OutVectorSInt("temnegnums",negnumbers,true);

	return 0;
}


int main__()
{
	_chdir("E:\\uiucCars\\CarData\\TrainImages");

	//PMStruc ptem;

	//ptem.loadFromAFile("firtem.txt");

	vector<featype> allfeas;
	allfeas.clear();
	allfeas=fileIOclass::InVector<featype>("allPostiveFeatures.txt");
	
	vector<vector<double> > TrmTx;
	TrmTx=fileIOclass::InVectorSDouble("../DimensionReduction.txt");

	TrmTx=keepFirSevDims(TrmTx,10);
	vector<vector<double > >  dats;
	dats.resize(allfeas.size(),vector<double>(0,0.0));
	for (int i=0;i<dats.size();i++)
	{
		dats[i]=allfeas[i].toVdouble();
	}
	
	ZeroMnVec(dats);
	dats=TransitMtx(dats,TrmTx);
	
	dats=addPositionsToData(dats,allfeas);
	//NormalVec(dats);
	fileIOclass::OutVectorSDouble("temp.txt",dats,true);

	PMStruc pedmd(5,PMStruc::normal,"pym");
	pedmd.generatePymFromdata(dats);
	pedmd.outToAFile("firtem.txt");

	
	return 0;

}


