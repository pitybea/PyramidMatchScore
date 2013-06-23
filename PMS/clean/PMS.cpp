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

PMStruc::PMStruc()
{
	totalLvls=-1;
	mymode=PMStruc::unset;
	//name="";
}
PMStruc::PMStruc(PyrMode p,string s)
{
	totalLvls=0;
	mymode=p;
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
		return MatchDttoPym(dataset,ExcluMode,scoreAllLevel);
		break;
	case PMStruc::average:
		return MatchDttoPymAv(dataset,ExcluMode,scoreAllLevel);
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
	case PMStruc::average:
		return dataToPymAver(data);
		break;
	default:
		return 0;
		break;
	}
}

void PMStruc::valueToInx(pair<double,double> minMax,pair<double,double>& aAndB,int levl)
{
	int totl=pow2[levl];
	minMax.first-=0.001;
	minMax.second+=0.001;
	aAndB.first=(double)(totl)/(minMax.second-minMax.first);
	aAndB.second=0-aAndB.first*minMax.first;
}

bool PMStruc::dataToPymLvl(vector<vector<double> > datas,int lvel,map<int,map<int,int> >& pymlvl,vector<pair<double,double> > aAndB)
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
			indtem=aAndB[j].first*datas[i][j]+aAndB[j].second;
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
			indtem=aAndB[j].first*datas[i][j]+aAndB[j].second;
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
	fprintf(fp,"%d %d\n%d\n",mymode,totalLvls,pym.size());
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
	case PMStruc::average:
		fprintf(fp,"%d\n",intvDecs.size());
		for (int i = 0; i < intvDecs.size(); i++)
		{
			printintvDecs(intvDecs[i],fp);
		}

		break;
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
	fscanf(fp,"%d %d\n%d\n",&mymode,&totalLvls,&temsz);
	
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
	case PMStruc::average:
		int ivDS;
		fscanf(fp,"%d\n",&ivDS);
		for (int i = 0; i < ivDS; i++)
		{
			vector<vector<double> > indsOne;
			indsOne=loadintvDecs(fp);
			intvDecs.push_back(indsOne);
		
		}

		break;
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
}

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

		/*
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
		}*/

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

		bool goGon=false;
		int goodi=0;
		while(!goGon&&goodi<6)
		{
			vector<pair<double,double> > aAb;
			aAb.resize(data[0].size(),pair<double,double>(0.0,0.0));

			for (int i=0;i<data[0].size();i++)
			{
				valueToInx(minmax[i],aAb[i],goodi);
			}
			map<int,map<int,int> > pymlv;
			goGon=dataToPymLvl(data,goodi,pymlv,aAb);
			pym.push_back(pymlv);

			aAbs.push_back(aAb);
			weights.push_back(pow2[goodi]);
			goodi+=1;
			totalLvls+=1;
		}
	
	
	}
	
	return 0;
}
int PMStruc:: matchDToOneLv(vector<vector<double> > dataset,int levl,map<int,map<int,int> > pmlv,vector<vector<double> > invs,bool ExcluMode )
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
	for (int i=0;i<dataset.size();i++)
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
					}
				}
				else
				{	
					if (pmlv[fisI][secI]>0)
					{
						res+=1;
						pmlv[fisI][secI]-=1;
					}
				}

			}
		}
	}
	return res;
}
int PMStruc:: matchDToOneLv(vector<vector<double> > dataset,int levl,map<int,map<int,int> > pmlv,vector<pair<double,double> > aAndB, bool ExcluMode )
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
	for (int i=0;i<dataset.size();i++)
	{
		int fisI(0),secI(0);
		for (int j=0;j<haldim;j++)
		{
			int indtem;
			indtem=aAndB[j].first*dataset[i][j]+aAndB[j].second;
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
			indtem=aAndB[j].first*dataset[i][j]+aAndB[j].second;
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
					}
				}
				else 
					{
						if (pmlv[fisI][secI]>0)
						{
							res+=1;
							pmlv[fisI][secI]-=1;
						}
				}
			}
		}
	}
	return res;
}

double PMStruc::MatchDttoPym(vector<vector<double> > dataset,bool ExcluMode,vector<int> & mnumbers)
{
	//vector<int> mnumbers;
	mnumbers.resize(pym.size(),0);
	for (int i=0;i<mnumbers.size();i++)
	{
		mnumbers[i]=matchDToOneLv(dataset,i,pym[i],aAbs[i],ExcluMode);
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
void setDoubleVec(vector<double>& to,vector<double> frm)
{
	to.clear();
	to.insert(to.end(),frm.begin(),frm.end());
}

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
}

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
}
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