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
PMStruc::PMStruc(int i)
{
	totalLvls=i;
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

void PMStruc::valueToInx(pair<double,double> minMax,pair<double,double>& aAndB,int levl)
{
	int totl=pow2[levl];
	minMax.first-=0.001;
	minMax.second+=0.001;
	aAndB.first=(double)(totl)/(minMax.second-minMax.first);
	aAndB.second=0-aAndB.first*minMax.first;
}

void PMStruc::dataToPymLvl(vector<vector<double> > datas,int lvel,map<int,map<int,int> >& pymlvl,vector<pair<double,double> > aAndB)
{
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
		if ((re>=inv[i])&&(re<inv[i+1]))
		{
			return re;
		}
	}
	return re;
}
void PMStruc::dataToPymLvl(vector<vector<double> > datas,int lvel,map<int,map<int,int> >& pymlvl,vector<vector<double> > aintvl)
{
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
		intvDecs.resize(totalLvls,vector<vector<double> >(dimension,tvd));
		
		initintvs(intvDecs,data);
		//for(int lvel=0;lvel<totalLvls;lvel++)
		//	for(int i=0;i<dimension;i++)
		//		dtinterv(data,intvDecs[lvel][i],i,lvel);

		map<int,map<int,int> > pymlv;
		pymlv.clear();
		pym.resize(totalLvls,pymlv);
		for (int i=0;i<totalLvls;i++)
		{
			dataToPymLvl(data,i,pym[i],intvDecs[i]);
		}

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
		aAbs.resize(totalLvls,vector<pair<double,double> >(data[0].size(),pair<double,double>(0.0,0.0)));
		for (int lvel=0;lvel<totalLvls;lvel++)
		{
			for (int i=0;i<data[0].size();i++)
			{
				//dtmMinMx(datas,i,minmax[i]);
				valueToInx(minmax[i],aAbs[lvel][i],lvel);
			}
		}
		map<int,map<int,int> > pymlv;
		pymlv.clear();
		pym.resize(totalLvls,pymlv);
		for (int i=0;i<totalLvls;i++)
		{
			dataToPymLvl(data,i,pym[i],aAbs[i]);
		}
	
	}
	
	return 0;
}
int PMStruc:: matchDToOneLv(vector<vector<double> > dataset,int levl,map<int,map<int,int> > pmlv,vector<vector<double> > invs )
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
				if (pmlv[fisI][secI]>0)
				{
					res+=1;
					pmlv[fisI][secI]=-1;
				}
			}
		}
	}
	return res;
}
int PMStruc:: matchDToOneLv(vector<vector<double> > dataset,int levl,map<int,map<int,int> > pmlv,vector<pair<double,double> > aAndB )
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
				if (pmlv[fisI][secI]>0)
				{
					res+=1;
					pmlv[fisI][secI]=-1;
				}
			}
		}
	}
	return res;
}

double PMStruc::MatchDttoPym(vector<vector<double> > dataset)
{
	vector<int> mnumbers;
	mnumbers.resize(pym.size(),0);
	for (int i=0;i<mnumbers.size();i++)
	{
		mnumbers[i]=matchDToOneLv(dataset,i,pym[i],aAbs[i]);
	}
	for (int i=0;i<mnumbers.size()-1;i++)
	{
		mnumbers[i]=mnumbers[i]-mnumbers[i+1];
	}
	vector<double> weights;
	weights.resize(pym.size(),0.0);
	double reslt(0.0);
	for (int i=0;i<pym.size();i++)
	{
		reslt+=mnumbers[i]*pow2[i]*dataset[0].size();
	}
	return reslt;
}

double PMStruc::MatchDttoPymAv(vector<vector<double> > dataset)
{
	vector<int> mnumbers;
	mnumbers.resize(pym.size(),0);
	for (int i=0;i<mnumbers.size();i++)
	{
		mnumbers[i]=matchDToOneLv(dataset,i,pym[i],intvDecs[i]);
	}
	for (int i=0;i<mnumbers.size()-1;i++)
	{
		mnumbers[i]=mnumbers[i]-mnumbers[i+1];
	}
	vector<double> weights;
	weights.resize(pym.size(),0.0);
	double reslt(0.0);
	for (int i=0;i<pym.size();i++)
	{
		reslt+=mnumbers[i]*pow2[i]*dataset[0].size();
	}
	return reslt;
}
