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



int main()
{
	_chdir("E:\\uiucCars\\CarData\\TrainImages");

	vector<featype> allfeas;
	allfeas.clear();
	allfeas=fileIOclass::InVector<featype>("allPostiveFeatures.txt");
	
	vector<vector<double> > TrmTx;
	TrmTx=fileIOclass::InVectorSDouble("../DimensionReduction.txt");

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

	PMStruc pedmd(5);
	pedmd.dataToPym(dats);

	return 0;

}


