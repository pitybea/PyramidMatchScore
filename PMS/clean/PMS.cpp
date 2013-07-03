#include "PMS.h"



static float pow2[]=
{
	1,
	2,
	4,
	8,
	16,
	32,
	64,
	128,
	256,
	512,
	1024,
	2048,
	4096,
	8192,
	16384,
	32768,
	65536,
	131072,
	262144,
	524288,
	1048576,
	2097152,
	4194304,
	8388608

};


PMStruc::PMStruc(PyrMode p)
{

	mymode=p;
	numOfData=0;

	
	//name=s;
	/*
	weights.resize(i,0.0);
	for(int i=0;i<weights.size();i++)
	{
		weights[i]=pow2[i];
	}*/
}

template<class T>
void dtmMinMx(vector<vector<T> > data,vector<pair<T,T> >& minmaxs)
{

	for (int ind=0;ind<data[0].size();ind++)
	{
		minmaxs[ind].first=data[0][ind];
		minmaxs[ind].second=data[0][ind];
		for (int i=0;i<data.size();i++)
		{
			if (data[i][ind]<minmaxs[ind].first)
			{
				minmaxs[ind].first=data[i][ind];
			}
			else if (data[i][ind]>minmaxs[ind].second)
			{
				minmaxs[ind].second=data[i][ind];
			}
		}
	}
	
}


double PMStruc::givePyramidMatchScore(vector<vector<double> > dataset,bool ExcluMode,vector<int> & scoreAllLevel)
{
	switch (mymode)
	{
	case PMStruc::normal:
		return MatchDttoPym(dataset,ExcluMode,scoreAllLevel,false);
		break;
/*	case PMStruc::average:
		return MatchDttoPymAv(dataset,ExcluMode,scoreAllLevel);
		break;*/
	case PMStruc::inverse:
		return MatchDttoPym(dataset,ExcluMode,scoreAllLevel,true);
		break;
	case PMStruc::postitionSpecific:
		return MatchPosDttoPym(dataset,ExcluMode,scoreAllLevel,false,posorder);
		break;
	default:
		return 0.0;
		break;
	}
}

int PMStruc::generatePymFromdata(vector<vector<double> > data)
{
	switch (mymode)
	{
	case PMStruc::normal:
		return dataToPym(data);
		break;
/*	case PMStruc::average:
		return dataToPymAver(data);
		break;*/
	case PMStruc::inverse:
		return dataToPym(data);
		break;
	case PMStruc::postitionSpecific:
		return dataToPym(data);
		break;
	default:
		return 0;
		break;
	}
}


pair<int,int> PMStruc::dataToTwoPosInx(int alvel,int plvel,vector<double> data)
{
	int dimension=data.size();
	int haldim=dimension/2;

	int fisI(0),secI(0);
	for (int j=0;j<haldim;j++)
		{
			int indtem;
			indtem=aAbs[alvel][j].first*data[j]+aAbs[alvel][j].second;
			if (indtem<0)
			{
				indtem=0;
			}
			if (indtem>(twExMs[alvel]-1))
			{
				indtem=twExMs[alvel]-1;
			}
			fisI+=twoExs[alvel][j]*indtem;
		}
		for (int j=haldim;j<dimension-2;j++)
		{
			int indtem;
			indtem=aAbs[alvel][j].first*data[j]+aAbs[alvel][j].second;
			if (indtem<0)
			{
				indtem=0;
			}
			if (indtem>(twExMs[alvel]-1))
			{
				indtem=twExMs[alvel]-1;
			}
			secI+=twoExs[alvel][j]*indtem;
		}
		for (int j = dimension-2; j < dimension; j++)
		{
			int indtem;
			indtem=aAbs[plvel][j].first*data[j]+aAbs[plvel][j].second;
			if (indtem<0)
			{
				indtem=0;
			}
			if (indtem>(twExMs[plvel]-1))
			{
				indtem=twExMs[plvel]-1;
			}
			secI+=twoExs[max(plvel,alvel)][j]*indtem;
		}

		return pair<int,int>(fisI,secI);
}


pair<int,int> PMStruc::dataToTwoIndx(int lvel,vector<double> data)
{

		int dimension=data.size();
		int haldim=dimension/2;

		int fisI(0),secI(0);
		for (int j=0;j<haldim;j++)
		{
			int indtem;
			indtem=aAbs[lvel][j].first*data[j]+aAbs[lvel][j].second;
			if (indtem<0)
			{
				indtem=0;
			}
			if (indtem>(twExMs[lvel]-1))
			{
				indtem=twExMs[lvel]-1;
			}
			fisI+=twoExs[lvel][j]*indtem;
		}
		for (int j=haldim;j<dimension;j++)
		{
			int indtem;
			indtem=aAbs[lvel][j].first*data[j]+aAbs[lvel][j].second;
			if (indtem<0)
			{
				indtem=0;
			}
			if (indtem>(twExMs[lvel]-1))
			{
				indtem=twExMs[lvel]-1;
			}
			secI+=twoExs[lvel][j]*indtem;
		}

		return pair<int,int> (fisI,secI);
}

void valueToInx(pair<double,double> minMax,pair<double,double>& aAndB,int levl)
{
	int totl=pow2[levl];
	minMax.first-=0.001;
	minMax.second+=0.001;
	aAndB.first=(double)(totl)/(minMax.second-minMax.first);
	aAndB.second=0-aAndB.first*minMax.first;
}



bool PMStruc::dataToPymLvl(vector<vector<double> > datas,int lvel,map<int,map<int,int> >& pymlvl)
{
	bool fineEngough=true;
	int dimension=datas[0].size();
	
	int haldim=datas[0].size()/2;


	for (int i=0;i<datas.size();i++)
	{
		pair<int,int> ss;
		if (mymode!=PMStruc::postitionSpecific)
			ss=dataToTwoIndx(lvel,datas[i]);
		else
			ss= dataToTwoPosInx(lvel%LevelLimit,lvel/LevelLimit,datas[i]);
		int fisI(ss.first),secI(ss.second);

		if (pymlvl.count(fisI)>0)
		{
			if (pymlvl[fisI].count(secI)>0)
			{
				pymlvl[fisI][secI]+=1;
				fineEngough=false;
			}
			else
			{
				pymlvl[fisI].insert(pair<int,int>(secI,1));
			}
		}
		else
		{
			map<int,int> temlv;
			temlv.insert(pair<int,int>(secI,1));
			pymlvl.insert(pair<int,map<int,int> >(fisI,temlv));
		}

	}
	return fineEngough;
}



int invdvalue(double a,vector<double> inv)
{
	int re=0;
	if (inv.size()==0)
	{
		return re;
	}
	if ((re<inv[0]))
	{
		return re;
	}
	if (inv.size()==1)
	{
		return (a>inv[0])?1:0;
	}
	for (int i=0;i<inv.size()-1;i++)
	{
		re+=1;
		if ((a>=inv[i])&&(a<inv[i+1]))
		{
			return re;
		}
	}
	return re;
}

map<int,map<int,int> > loadPyramidLv(FILE* fp)
{
	map<int,map<int,int> > result;
	int pymlvS;
	fscanf(fp,"%d\n",&pymlvS);
	for (int i = 0; i < pymlvS; i++)
	{
		int a,bS;
		fscanf(fp,"%d %d\n",&a,&bS);
		map<int,int> b;
		for (int j = 0; j < bS; j++)
		{
			int ta,tb;
			fscanf(fp,"%d %d\t",&ta,&tb);
			b.insert(pair<int,int>(ta,tb));
		}

		result.insert(pair<int,map<int,int> >(a,b));
		fscanf(fp,"\n");
	}
	return result;
}

void printPyramidLv(map<int,map<int,int> > pymlv,FILE* fp)
{
	fprintf(fp,"%d\n",pymlv.size());
	for(map<int,map<int,int> >::iterator ii=pymlv.begin(); ii!=pymlv.end(); ++ii)
	{
		int a=(*ii).first;
		map<int,int> b=(*ii).second;
		fprintf(fp,"%d %d\n",a,b.size());
		for(map<int,int>::iterator ij=b.begin();ij!=b.end();++ij)
		{
			fprintf(fp,"%d %d\t",(*ij).first,(*ij).second);
		}
		fprintf(fp,"\n");
	}

}

vector<pair<double,double> > loadAandBsOne(FILE* fp)
{
	vector<pair<double,double> > result;
	int abS;
	fscanf(fp,"%d\n",&abS);
	for (int i = 0; i < abS; i++)
	{
		double ta,tb;
		fscanf(fp,"%lf %lf\t",&ta,&tb);
		result.push_back(pair<double,double>(ta,tb));
	}
	return result;
}

void printAandBs(vector<pair<double,double> > abs,FILE* fp)
{
	fprintf(fp,"%d\n",abs.size());
	for(int i=0;i<abs.size();i++)
	{
		fprintf(fp,"%lf %lf\t",abs[i].first,abs[i].second);
	}
}

vector<vector<double> > loadintvDecs(FILE* fp)
{
	vector<vector<double> > result;
	int rSz;
	fscanf(fp,"%d\n",&rSz);
	for (int i = 0; i < rSz; i++)
	{
		int sSz;
		fscanf(fp,"%d\n",&sSz);
		vector<double> tvd;
		tvd.clear();
		tvd.resize(sSz,0.0);
		for (int j = 0; j < sSz; j++)
		{
//			double td;
			fprintf(fp,"%lf\t",tvd[j]);
		}
		fscanf(fp,"\n");
		result.push_back(tvd);
	}
	return result;
}
void printintvDecs(vector<vector<double> > intvecs,FILE* fp)
{
	fprintf(fp,"%d\n",intvecs.size());
	for(int i=0;i<intvecs.size();i++)
	{
		fprintf(fp,"%d\n",intvecs[i].size());
		for (int j = 0; j < intvecs[i].size(); j++)
		{
			fprintf(fp,"%lf\t",intvecs[i][j]);
		}
		fprintf(fp,"\n");
	}

}
void PMStruc::outToAFile(string filename)
{
	FILE* fp;
	fp=fopen(filename.c_str(),"w");
	fprintf(fp,"%d\n%d\n",mymode,pym.size());
	for(int i=0;i<pym.size();i++)
	{
		printPyramidLv(pym[i],fp);
		
	}

	switch (mymode)
	{
	case PMStruc::normal:
		fprintf(fp,"%d\n",aAbs.size());
		for (int i = 0; i < aAbs.size(); i++)
		{
			printAandBs(aAbs[i],fp);
			fprintf(fp,"\n");
		}
		break;
/*	case PMStruc::average:
		fprintf(fp,"%d\n",intvDecs.size());
		for (int i = 0; i < intvDecs.size(); i++)
		{
			printintvDecs(intvDecs[i],fp);
		}

		break;*/
	default:
		break;
	}

	


	
	fclose(fp);



	FILE* ft;
	string wFname=filename+"_wgt.txt";
	ft=fopen(wFname.c_str(),"w");
	fprintf(ft,"%d\n",weights.size());
	for(int i=0;i<weights.size();i++)
	{
		fprintf(ft,"%lf ",weights[i]);
	}
	fclose(ft);
}

void PMStruc::loadFromAFile(string filename)
{
	FILE* fp;
	fp=fopen(filename.c_str(),"r");
	char tems[100];
	int temsz;
	fscanf(fp,"%d\n%d\n",&mymode,&temsz);
	
	pym.clear();
	for (int i = 0; i < temsz; i++)
	{
		map<int,map<int,int> > pymlv=loadPyramidLv(fp);
		
		pym.push_back(pymlv);
	}


	switch (mymode)
	{
	case PMStruc::normal:
		int abbS;
		fscanf(fp,"%d\n",&abbS);
		for (int i = 0; i < abbS; i++)
		{
			vector<pair<double,double> > abbOne;
			abbOne=loadAandBsOne(fp);
			aAbs.push_back(abbOne);
		//	printAandBs(aAbs[i],fp);
			fscanf(fp,"\n");
		}
		break;
	/*case PMStruc::average:
		int ivDS;
		fscanf(fp,"%d\n",&ivDS);
		for (int i = 0; i < ivDS; i++)
		{
			vector<vector<double> > indsOne;
			indsOne=loadintvDecs(fp);
			intvDecs.push_back(indsOne);
		
		}

		break;*/
	default:
		break;
	}

	fclose(fp);
	FILE* ft;
	string wFname=filename+"_wgt.txt";
	ft=fopen(wFname.c_str(),"r");
	int wSiz;
	fscanf(ft,"%d\n",&wSiz);
	weights.resize(wSiz,0.0);
	for(int i=0;i<wSiz;i++)
	{
		fscanf(ft,"%lf ",&weights[i]);
	}
	fclose(ft);


}
/*
bool PMStruc::dataToPymLvl(vector<vector<double> > datas,int lvel,map<int,map<int,int> >& pymlvl,vector<vector<double> > aintvl)
{
	bool fineEngough=true;

	int dimension=datas[0].size();
	int haldim=datas[0].size()/2;
	int twExM=pow2[lvel];

	vector<int> twoEx;
	twoEx.resize(dimension,0);
	for (int i=0;i<haldim;i++)
	{
		twoEx[i]=twoEx[i+haldim]=pow(pow2[lvel],i);
	}
	for (int i=0;i<datas.size();i++)
	{
		int fisI(0),secI(0);
		for (int j=0;j<haldim;j++)
		{
			int indtem;

			indtem=invdvalue(datas[i][j],aintvl[j]);
			
			if (indtem<0)
			{
				indtem=0;
			}
			if (indtem>(twExM-1))
			{
				indtem=twExM-1;
			}
			fisI+=twoEx[j]*indtem;
		}
		for (int j=haldim;j<dimension;j++)
		{
			int indtem;
			indtem=invdvalue(datas[i][j],aintvl[j]);
			if (indtem<0)
			{
				indtem=0;
			}
			if (indtem>(twExM-1))
			{
				indtem=twExM-1;
			}
			secI+=twoEx[j]*indtem;
		}
		if (pymlvl.count(fisI)>0)
		{
			if (pymlvl[fisI].count(secI)>0)
			{
				pymlvl[fisI][secI]+=1;
				fineEngough=false;
				//	cout<<pymlvl[fisI][secI]<<endl;
			}
			else
			{
				pymlvl[fisI].insert(pair<int,int>(secI,1));
			}
		}
		else
		{
			map<int,int> temlv;
			temlv.insert(pair<int,int>(secI,1));
			pymlvl.insert(pair<int,map<int,int> >(fisI,temlv));
		}

	}
	return fineEngough;
}*/

void dtinterv(vector<vector<double> > dataset,vector<double> & intvDec,int i,int goodi)
{
		vector<double> temvec;
		temvec.resize(dataset.size(),0.0);
		vector<int> inde_;
		inde_.clear();
		inde_.resize(dataset.size(),0);
		for (int j=0;j<dataset.size();j++)
		{
			temvec[j]=dataset[j][i];
			inde_[j]=j;
		}
		prshl(temvec,temvec.size(),inde_);

		int jj=goodi;
	
		int tot=pow2[jj];
		intvDec.clear();
		if (tot>1)
		{
			int step=dataset.size()/tot;
			for (int k=1;k<pow2[jj];k++)
			{
				intvDec.push_back(temvec[inde_[k*step]]);
			}
							
		}
		
}

int initintvs(vector<vector<vector<double> > > &intvDecs,vector<vector<double> > dataset )
{
	int dimension=dataset[0].size();
	for (int i=0;i<dimension;i++)
	{
		vector<double> temvec;
		temvec.resize(dataset.size(),0.0);
		vector<int> inde_;
		inde_.clear();
		inde_.resize(dataset.size(),0);
		for (int j=0;j<dataset.size();j++)
		{
			temvec[j]=dataset[j][i];
			inde_[j]=j;
		}
		prshl(temvec,temvec.size(),inde_);
		for (int j=0;j<intvDecs.size();j++)
		{
			int tot=pow2[j];
			intvDecs[j][i].clear();
			if (tot>1)
			{
				int step=dataset.size()/tot;
				for (int k=1;k<pow2[j];k++)
				{
					intvDecs[j][i].push_back(temvec[inde_[k*step]]);
				}
							
			}
		}
	}
	return 0;
}

/*
int PMStruc::dataToPymAver(vector<vector<double> > data)
{
	if (data.size()>0)
	{
		int dimension=data[0].size();
		vector<double> tvd;
		tvd.clear();

		vector<vector<vector<double> > > intvDecs_;
		intvDecs_.resize(10,vector<vector<double> >(dimension,tvd));
		
		initintvs(intvDecs_,data);
		bool goGon=false;
		int goodi=0;
		while(!goGon&&goodi<10)
		{
			
				//dtinterv(vector<vector<double> > data,vector<vector<double> >& intvDec,int i,int goodi);
			intvDecs.push_back(intvDecs_[goodi]);

			
			map<int,map<int,int> > pymlv;
			pymlv.clear();
			
			goGon=dataToPymLvl(data,goodi,pymlv,intvDecs_[goodi]);
			pym.push_back(pymlv);

			weights.push_back(pow2[goodi]);
			goodi+=1;
			totalLvls+=1;
		}

		//comment
		intvDecs.resize(totalLvls,vector<vector<double> >(dimension,tvd));
		
		initintvs(intvDecs,data);
		for(int lvel=0;lvel<totalLvls;lvel++)
			for(int i=0;i<dimension;i++)
				dtinterv(data,intvDecs[lvel][i],i,lvel);

		map<int,map<int,int> > pymlv;
		pymlv.clear();
		pym.resize(totalLvls,pymlv);
		for (int i=0;i<totalLvls;i++)
		{
			dataToPymLvl(data,i,pym[i],intvDecs[i]);
		}//comment
		numOfData+=data.size();
	}

	return 0;
}
*/


int PMStruc::GeneratePosWeightWithParameter(double ratio)
{
	for (int i = 0; i < LevelLimit*LevelLimit; i++)
	{
		int pi(i/LevelLimit);int ai(i%LevelLimit);
		weights[i]= pow((double)pow2[pi],ratio)*pow( pow2[ai],1-ratio);

	}

	return 0;
}

int PMStruc::dataToPym(vector<vector<double> > data)
{
	if (data.size()>0)
	{
		int dimension=data[0].size();
		vector<pair<double,double> > minmax;
		minmax.resize(data[0].size(),pair<double,double>(0.0,0.0));
		dtmMinMx(data,minmax);

	
		int haldim=dimension/2;
		
		for (int goodi = 0; goodi < LevelLimit; goodi++)
		{
			int twExM=pow2[goodi];
			vector<int> twoEx;
			twoEx.resize(dimension,0);
			for (int ti=0;ti<haldim;ti++)
			{
				twoEx[ti]=twoEx[ti+haldim]=pow(pow2[goodi],ti);
			}
			
			twExMs.push_back(twExM);

			twoExs.push_back(twoEx);

			vector<pair<double,double> > aAb;
			aAb.resize(data[0].size(),pair<double,double>(0.0,0.0));

			for (int i=0;i<data[0].size();i++)
			{
				valueToInx(minmax[i],aAb[i],goodi);
			}

			aAbs.push_back(aAb);
			
	

		}

		if(mymode!=PMStruc::postitionSpecific)
			for (int goodi = 0; goodi < LevelLimit; goodi++)
			{
				weights.push_back(pow2[goodi]);
				map<int,map<int,int> > pymlv;
				dataToPymLvl(data,goodi,pymlv);
				pym.push_back(pymlv);
			}
		else
		{
			weights.resize(LevelLimit*LevelLimit,0.0);
		
			for (int i = 0; i < LevelLimit*LevelLimit; i++)
			{

					int ai(i/LevelLimit);int aj(i%LevelLimit);
					map<int,map<int,int> > pymlv;
					dataToPymLvl(data,i,pymlv);
					pym.push_back(pymlv);
					weights[i]=sqrt((double)(pow2[ai]*pow2[aj]));


			}

		}

		numOfData+=data.size();
	
	}
	
	return 0;
}

/*
int PMStruc:: matchDToOneLv(vector<vector<double> > &dataset,int levl,map<int,map<int,int> > pmlv,vector<vector<double> > invs,bool ExcluMode )
{
	int haldim=dataset[0].size()/2;
	int dimension=dataset[0].size();
	int twExM=pow2[levl];

	vector<int> twoEx;
	twoEx.resize(dimension,0);
	for (int i=0;i<haldim;i++)
	{
		twoEx[i]=twoEx[i+haldim]=pow(pow2[levl],i);
	}

	int res(0);

	int siz=dataset.size();
	for (int i=siz-1;i>=0;i--)
	{
		int fisI(0),secI(0);
		for (int j=0;j<haldim;j++)
		{
			int indtem;
			//indtem=aAndB[j].first*dataset[i][j]+aAndB[j].second;
			indtem= invdvalue(dataset[i][j],invs[j]);
			if (indtem<0)
			{
				indtem=0;
			}
			if (indtem>(twExM-1))
			{
				indtem=twExM-1;
			}
			fisI+=twoEx[j]*indtem;
		}
		for (int j=haldim;j<dimension;j++)
		{
			int indtem;
			indtem= invdvalue(dataset[i][j],invs[j]);
			if (indtem<0)
			{
				indtem=0;
			}
			if (indtem>(twExM-1))
			{
				indtem=twExM-1;
			}
			secI+=twoEx[j]*indtem;
		}


		if (pmlv.count(fisI)>0)
		{
			if (pmlv[fisI].count(secI)>0)
			{
				if(ExcluMode)
				{
					if (pmlv[fisI][secI]>1)
					{
						res+=1;
						pmlv[fisI][secI]-=1;
						dataset.erase(dataset.begin()+i);
					}
				}
				else
				{	
					if (pmlv[fisI][secI]>0)
					{
						res+=1;
						pmlv[fisI][secI]-=1;
						dataset.erase(dataset.begin()+i);
					}
				}

			}
		}
	}
	return res;
}

*/
int PMStruc:: matchDToOneLv(vector<vector<double> >& dataset,int levl,map<int,map<int,int> > pmlv,bool ExcluMode )
{
	int haldim=dataset[0].size()/2;
	int dimension=dataset[0].size();

	int siz=dataset.size();
	int res(0);
	for (int i=siz-1;i>=0;i--)
	{
		
		pair<int,int> ss;
		if(mymode!=PMStruc::postitionSpecific)
			ss=dataToTwoIndx(levl,dataset[i]);
		else
			ss=dataToTwoPosInx(levl%LevelLimit,levl/LevelLimit,dataset[i]);

		int fisI(ss.first),secI(ss.second);

		if (pmlv.count(fisI)>0)
		{
			if (pmlv[fisI].count(secI)>0)
			{
				if(ExcluMode)
				{
					if (pmlv[fisI][secI]>1)
					{
						res+=1;
						pmlv[fisI][secI]-=1;
						dataset.erase(dataset.begin()+i);
					}
				}
				else 
					{
						if (pmlv[fisI][secI]>0)
						{
							res+=1;
							pmlv[fisI][secI]-=1;
							dataset.erase(dataset.begin()+i);
						}
				}
			}
		}
	}
	return res;
}




double PMStruc::MatchDttoPym(vector<vector<double> > dataset,bool ExcluMode,vector<int> & mnumbers,bool inverse)
{
	//vector<int> mnumbers;

	int siz=dataset.size();

	int dim=dataset[0].size();
	mnumbers.resize(pym.size(),0);
	for (int i=mnumbers.size()-1;i>=0;i--)
	{
		mnumbers[i]=matchDToOneLv(dataset,i,pym[i],ExcluMode);
	}
	
	//vector<double> weights;
	//weights.resize(pym.size(),0.0);
	double reslt(0.0);
	for (int i=0;i<pym.size();i++)
	{
		if(!inverse)
			reslt+=mnumbers[i]*weights[i]*dim;
		else
			reslt+=mnumbers[i]*(1.0/weights[i])*dim;
	}
	reslt/=siz;
	if(inverse)
	{
		
		reslt=1.0/reslt;
	}

		
	return reslt;
}

double PMStruc::MatchPosDttoPym(vector<vector<double> > dataset,bool ExcluMode,vector<int> & mnumbers,bool inverse,int order[LevelLimit*LevelLimit][2])
{
	//vector<int> mnumbers;
	mnumbers.resize(LevelLimit*LevelLimit,0);
	int dim=dataset[0].size();
	int siz=dataset.size();

	double reslt(0.0);
	for (int wi = 0; wi < LevelLimit*LevelLimit; wi++)
	{
		int pi=order[wi][1];
		int ai=order[wi][0];
		mnumbers[wi]=matchDToOneLv(dataset,pi*LevelLimit+ai,pym[pi*LevelLimit+ai],ExcluMode);
		if(!inverse)
			reslt+=mnumbers[wi]*weights[pi*LevelLimit+ ai]*dim;
		else
			reslt+=mnumbers[wi]*(1.0/weights[pi*LevelLimit+ ai ])*dim;
			

	}
		
	
	
	reslt/=siz;
	if(inverse)
	{
		
		reslt=1.0/reslt;
	}
		
	return reslt;
}




/*

double PMStruc::MatchDttoPymAv(vector<vector<double> > dataset,bool ExcluMode,vector<int> & mnumbers)
{
	//vector<int> mnumbers;
	mnumbers.resize(pym.size(),0);
	for (int i=0;i<mnumbers.size();i++)
	{
		mnumbers[i]=matchDToOneLv(dataset,i,pym[i],intvDecs[i],ExcluMode);
	}
	for (int i=0;i<mnumbers.size()-1;i++)
	{
		mnumbers[i]=mnumbers[i]-mnumbers[i+1];
	}
	//vector<double> weights;
	//weights.resize(pym.size(),0.0);
	double reslt(0.0);
	for (int i=0;i<pym.size();i++)
	{
		reslt+=mnumbers[i]*weights[i]*dataset[0].size();
	}
	return reslt;
}
*/

double multiIVecDvec(vector<int> ivec,vector<double> dvec)
{
	assert(ivec.size()==dvec.size());
	double result=0.0;
	for (int i = 0; i < ivec.size(); i++)
	{
		result+=ivec[i]*dvec[i];
	}
	return result;
}
/*
int PMStrainer::DecideBestScore(vector<double> temWeight,vector<vector<int> > posnums,vector<vector<int> > negnums,double& divValue)
{
	if(posnums.size()==0||negnums.size()==0)
		return -1;

	assert(temWeight.size()==posnums[0].size());
	assert(temWeight.size()==negnums[0].size());

	vector<double> allWeights;
	allWeights.resize(posnums.size()+negnums.size(),0.0);
	vector<int> orderIdx;
	orderIdx.resize(posnums.size()+negnums.size(),0);

	for (int i = 0; i < posnums.size(); i++)
	{
		allWeights[i]=multiIVecDvec(posnums[i],temWeight);
		orderIdx[i]=i;
	}
	int tn=posnums.size();
	for (int i = 0; i < negnums.size(); i++)
	{
		allWeights[tn+i]=multiIVecDvec(negnums[i],temWeight);
		orderIdx[tn+i]=tn+i;
	}
	
	vector<int> allLabels;
	vector<int> tposL;
	tposL.resize(posnums.size(),1);
	allLabels.insert(allLabels.end(),tposL.begin(),tposL.end());
	
	vector<int> nposL;
	nposL.resize(negnums.size(),-1);
	allLabels.insert(allLabels.end(),nposL.begin(),nposL.end());

	prshl(allWeights,allWeights.size(),orderIdx);
	vector<int> PosCul;
	vector<int> NegCul;
	PosCul.resize(orderIdx.size(),0);
	NegCul.resize(orderIdx.size(),0);

	for(int i=0;i<orderIdx.size();i++)
	{
		int tlabel=allLabels[orderIdx[i]];
		if (i==0)
		{
			if(tlabel>0)
			{
				PosCul[0]=1;
				NegCul[0]=0;
			}
			else
			{
				PosCul[0]=0;
				NegCul[0]=1;
			}
		}
		else
		{
			if(tlabel>0)
			{
				PosCul[i]=PosCul[i-1]+1;
				NegCul[i]=NegCul[i-1];
			}
			else
			{
				PosCul[i]=PosCul[i-1];
				NegCul[i]=NegCul[i-1]+1;
			}
		}
	}
	int allposN,allnegN;
	allposN=posnums.size();
	allnegN=negnums.size();

	int bestPosition=-1;
	int bestNum=-1;

	for (int i = 0; i < orderIdx.size()-1; i++)
	{
		int corNum=NegCul[i]+allposN-PosCul[i]+1;
		if(corNum>bestNum)
		{
			bestPosition=i;
			bestNum=corNum;
		}
	}
	divValue=(allWeights[orderIdx[bestPosition]]+allWeights[orderIdx[bestPosition+1]])/2;
	return bestNum;
}
*/

/*	
PMStrainer::PMStrainer(PMStruc pm,TrainStrategy s,double d,double l,string ts)
{
	myStrgy=s;
	pmStrct=pm;
	stepWidth=d;
	largeBand=l;
	name=ts;
}

void PMStrainer::Train(vector<vector<int> > posnums,vector<vector<int> > negnums)
{
	


	switch (myStrgy)
	{
	case PMStrainer::allAtOnce:
		trainAllAtOnce(posnums,negnums);
		break;
	case PMStrainer::oneAtOnce:
		trainOneAtOnce(posnums,negnums);
		break;
	default:
		break;
	}
}
*/
void setDoubleVec(vector<double>& to,vector<double> frm)
{
	to.clear();
	to.insert(to.end(),frm.begin(),frm.end());
}
/*
void PMStrainer::trainOneAtOnce(vector<vector<int> > posnums,vector<vector<int> > negnums)
{
	vector<double> weights;
	weights.resize(pmStrct.weights.size(),0.0);
	weights[0]=1.0;
	int bNum(-1);
	double bScore(-1.0);
	vector<double> bWeights;
	bWeights.clear();

	int tnum;
	double tvalue;
	for (int i = 1; i < weights.size(); i++)
	{
		double prW=weights[i-1];
		double hBd=pmStrct.weights[i]*largeBand;
		for(double tw=prW;tw< hBd;tw+=stepWidth)
		{
			tnum=DecideBestScore(weights,posnums,negnums,tvalue);
			if(tnum>bNum)
			{
				bNum=tnum;
				bScore=tvalue;
				setDoubleVec(bWeights,weights);
			}
		}
	}

	bestNumber=bNum;
	bestDivScore=bScore;
	setDoubleVec(bestWeights,bWeights);
}*/

bool endWhile(vector<double> jg,vector<double> tojg)
{
	assert(jg.size()==tojg.size());
	for (int i = 0; i < jg.size(); i++)
	{
		if(tojg[i]<jg[i])
			return false;
	}
	return true;
}


void nextweights(vector<double> bd,vector<double>& wts,double stw)
{
	int idx=wts.size()-1;
	bool stpPr=false;
	while(!stpPr)
	{
		wts[idx]+=stw;
		if(wts[idx]<bd[idx])
		{
			stpPr=true;
			if(idx!=(wts.size()-1))
			{
				for(int i=idx+1;i<wts.size();i++)
					wts[i]=wts[idx];							
			}
		}
		else
		{
			idx-=1;
		}
		if(idx==0)
		{
			stpPr=true;
			break;
		}

		
	}
}
/*
void PMStrainer::trainAllAtOnce(vector<vector<int> > posnums,vector<vector<int> > negnums)
{

	vector<double> weights;
	weights.resize(pmStrct.weights.size(),1.0);
	//weights[0]=1.0;
	int bNum(-1);
	double bScore(-1.0);
	vector<double> bWeights;
	bWeights.clear();

	int tnum;
	double tvalue;

	vector<double> allBds;
	allBds.resize(weights.size(),0.0);
	allBds[0]=0.9;

	for (int i = 1; i < weights.size(); i++)
	{
		allBds[i]=pmStrct.weights[i]*largeBand;
	}
	while(!endWhile(allBds,weights))
	{
		
		tnum=DecideBestScore(weights,posnums,negnums,tvalue);
		if(tnum>bNum)
		{
			bNum=tnum;
			bScore=tvalue;
			setDoubleVec(bWeights,weights);
		}
		nextweights(allBds,weights,stepWidth);
	}

	bestNumber=bNum;
	bestDivScore=bScore;
	setDoubleVec(bestWeights,bWeights);
}

void PMStrainer::OutTofile()
{
	string s=name+"_bestWts.txt";
	FILE* fp;
	fp=fopen(s.c_str(),"w");
	fprintf(fp,"%d %lf\n%d\n",bestNumber,bestDivScore,bestWeights.size());
	for (int i = 0; i < bestWeights.size(); i++)
	{
		fprintf(fp,"%lf\n",bestWeights[i]);
	}
	fclose(fp);
}*/
/*
double prW=weights[i-1];
		double hBd=pmStrct.weights[i]*largeBand;
		for(double tw=prW;tw< hBd;tw+=stepWidth)
		{
			tnum=DecideBestScore(weights,posnums,negnums,tvalue);
			if(tnum>bNum)
			{
				bNum=tnum;
				bScore=tvalue;
				setDoubleVec(bWeights,weights);
			}
		}

*/

PMSEnsemble::PMSEnsemble()
{
	pyms.clear();
}

int PMSEnsemble::generateAaBsFromdata(vector<vector<double> > data)
{
	if (data.size()>0)
	{
		dimension=data[0].size();
		vector<pair<double,double> > minmax;
		minmax.resize(data[0].size(),pair<double,double>(0.0,0.0));
		dtmMinMx(data,minmax);
	
		int goodi=0;
		while(goodi<LevelLimit)
		{
			vector<pair<double,double> > aAb;
			aAb.resize(data[0].size(),pair<double,double>(0.0,0.0));
			for (int i=0;i<data[0].size();i++)
			{
				valueToInx(minmax[i],aAb[i],goodi);
			}			
			aAbs.push_back(aAb);	
			goodi+=1;

		}
	
	
	}
	
	return 0;
};

int PMStruc::initPymWithABs(vector<vector<pair<double,double> > > abS,int dimension)
{
	if (abS.size()>0)
	{


	
		int haldim=dimension/2;

		aAbs.clear();
		aAbs.insert(aAbs.end(),abS.begin(),abS.end());

		twoExs.clear();
		twExMs.clear();

		for (int i = 0; i < aAbs.size(); i++)
		{
			map<int,map<int,int> > pymlv;	
			pym.push_back(pymlv);
			weights.push_back(pow2[i]);
//			totalLvls+=1;



			
			int twExM=pow2[i];
			vector<int> twoEx;
			twoEx.resize(dimension,0);
			for (int ti=0;ti<haldim;ti++)
			{
				twoEx[ti]=twoEx[ti+haldim]=pow(pow2[i],ti);
			}
			twExMs.push_back(twExM);
			twoExs.push_back(twoEx);
		}

		

	
	
	}
	
	return 0;
};
int PMStruc::AddoneData(vector<double> data,bool AddOrMinus)
{

	
	int dimension=data.size();
	int haldim=dimension/2;

	for (int levi = 0; levi < pym.size(); levi++)
	{
		int lvel=levi;
	
		int fisI(0),secI(0);
		for (int j=0;j<haldim;j++)
		{
			int indtem;
			indtem=aAbs[levi][j].first*data[j]+aAbs[levi][j].second;
			if (indtem<0)
			{
				indtem=0;
			}
			if (indtem>(twExMs[levi]-1))
			{
				indtem=twExMs[levi]-1;
			}
			fisI+=twoExs[levi][j]*indtem;
		}
		for (int j=haldim;j<dimension;j++)
		{
			int indtem;
			indtem=aAbs[levi][j].first*data[j]+aAbs[levi][j].second;
			if (indtem<0)
			{
				indtem=0;
			}
			if (indtem>(twExMs[levi]-1))
			{
				indtem=twExMs[levi]-1;
			}
			secI+=twoExs[levi][j]*indtem;
		}

		
		if (pym[levi].count(fisI)>0)
		{
			if (pym[levi][fisI].count(secI)>0)
			{
				if(AddOrMinus)
					pym[levi][fisI][secI]+=1;
				else
					pym[levi][fisI][secI]-=1;
				//fineEngough=false;
				//	cout<<pymlvl[fisI][secI]<<endl;
			}
			else
			{
				if(AddOrMinus)
					pym[levi][fisI].insert(pair<int,int>(secI,1));
				else
					pym[levi][fisI].insert(pair<int,int>(secI,-1));
			}
		}
		else
		{
			map<int,int> temlv;
			if(AddOrMinus)
				temlv.insert(pair<int,int>(secI,1));
			else
				temlv.insert(pair<int,int>(secI,-1));
			pym[levi].insert(pair<int,map<int,int> >(fisI,temlv));
		}
	}
	

	numOfData+=1;
	return 0;
};



int PMStruc::AddSeverlData(vector<vector<double> > data,bool AddorMinus)
{
	if(data.size()>0)
		for(auto d :data)
			AddoneData(d,AddorMinus);

	return 0;
}


int PMSEnsemble::generateStructureFromData(vector<vector<vector<double> > > data)
{

	if (data.size()>0)
		for(auto d :data)
			if(pyms.size()==0)
			{
				PMStruc tp(PMStruc::normal);
				tp.initPymWithABs(aAbs,dimension);
				tp.AddSeverlData(d,true);
				pyms.push_back(tp);
			}
			else
			{
				vector<double> tr;
				tr.resize(pyms.size(),0.0);
				for(int i=0;i<tr.size();i++)
				{
					vector<int> temiv;
					tr[i]=pyms[i].givePyramidMatchScore(d,false,temiv)/sqrt((double)pyms[i].getNumofData());
					
				}
				auto sd=maxminValAndInx(tr,true);
				if(sd.second>threshold)
					pyms[sd.first].AddSeverlData(d,true);
				else
				{
					
					PMStruc tp(PMStruc::normal);
					tp.initPymWithABs(aAbs,dimension);
					tp.AddSeverlData(d,true);
					pyms.push_back(tp);
				}


			}

	

	return 0;
}

double PMSEnsemble::givePyramidMatchScore(vector<vector<double> > dataset)
{
	if(pyms.size()==0)
		return 0.0;
	else
	{
		vector<double> rslts;
		rslts.clear();
		for(auto s:pyms)
			rslts.push_back(s.givePyramidMatchScore(dataset,false,vector<int>(0,0)));

		return maxminValAndInx(rslts,true).second;
	}
}


trainVecs:: trainVecs(trainVecs::strategies str,vector<vector<double> > lhb,vector<double> st,vector<vector<double> > pv,vector<vector<double> > nv)
{
	mystr=str;
	lowHighBands=lhb;
	stps=st;
	posvecs=pv;
	negvecs=nv;
}

void trainVecs::train()
{
	switch (mystr)
	{
	case trainVecs::allAtOnce:
		trainAllAtOnce();
		break;
	case trainVecs::oneAtonce:
		trainOneAtOnce();
		break;
	default:
		break;
	}
}

void trainVecs::trainOneAtOnce()
{
}

vector<double> nextweights(vector<vector <double> > bds, vector<double> owts,vector<double> stws,int &which)
{

	vector<double> wts(owts);

	int idx=wts.size()-1;
	bool stpPr=false;
	while(!stpPr)
	{
		wts[idx]+=stws[idx];
		if(wts[idx]<bds[1][idx])
		{
			stpPr=true;
			which=idx;
			if(idx!=(wts.size()-1))
			{
				for(int i=idx+1;i<wts.size();i++)
					wts[i]=bds[0][i];							
			}

		}
		else
		{
			idx-=1;
		}
		if(idx==0)
		{
			stpPr=true;
			break;
		}

		
	}
	return wts;
}


int DecideBestScore(vector<double> allWeights,vector<int> lblIdx,int allp,int alln)
{
	
	prshl(allWeights,allWeights.size(),lblIdx);
	int bestNum=-1;
	int nnum=0;
	int pnum=0;
	
	for (int i = 0; i < lblIdx.size(); i++)
	{
		if (lblIdx[i]==1)
		{
			pnum+=1;
		}
		else
		{
			nnum+=1;
		}
		
		int crrctNum=nnum+allp-pnum;

		if(crrctNum>bestNum)
		{
			bestNum=crrctNum;
		}


	}

	return bestNum;
}

double multiIVecDvec(vector<double> ivec,vector<double> dvec,int frm)
{

	double result=0.0;
	for (int i = 0; i < dvec.size(); i++)
	{
		result+=ivec[frm+i]*dvec[i];
	}
	return result;
}

void trainVecs::trainAllAtOnce()
{
	vector<double> tweights;
	tweights.clear();
	tweights.insert(tweights.end(),lowHighBands[0].begin(),lowHighBands[0].end());

	vector<double> allScores;
	allScores.resize(posvecs.size()+negvecs.size());

	int allp=posvecs.size();
	int alln=negvecs.size();

	

	vector<vector<double> > allvects;
	allvects.clear();
	allvects.insert(allvects.end(),posvecs.begin(),posvecs.end());
	allvects.insert(allvects.end(),negvecs.begin(),negvecs.end());

	vector<int> lblIdx;
	lblIdx.clear();
	vector<int> plbs(allp,1);
	vector<int> nlbs(alln,-1);
	lblIdx.insert(lblIdx.end(),plbs.begin(),plbs.end());
	lblIdx.insert(lblIdx.end(),nlbs.begin(),nlbs.end());

	for (int i = 0; i < allvects.size(); i++)
	{
		allScores[i]=multiIVecDvec(allvects[i],tweights,0);
		
	}
	
	vector<double> newweights;

	int dim=posvecs[0].size();

	int bestnum=-1;

	while(!endWhile(lowHighBands[1],tweights))
	{
		int inx;
		newweights=nextweights(lowHighBands,tweights,stps,inx);
		
		vector<double> diff;
		diff.resize(dim-inx,0.0);
		
		for (int i = 0; i < diff.size(); i++)
		{
			diff[i]=newweights[inx+i]-tweights[inx+i];
		}
		
		if (inx<4)
		{
			printf("%d ",inx);
		}

		for (int i = 0; i < allScores.size(); i++)
		{
			 allScores[i]+=multiIVecDvec(allvects[i],diff,inx);
		}

		int tnum=DecideBestScore(allScores,lblIdx,allp,alln);
		
		if (tnum>bestnum)
		{
			setDoubleVec(bstWeights,newweights);
		}

		setDoubleVec(tweights,newweights);

	}

}