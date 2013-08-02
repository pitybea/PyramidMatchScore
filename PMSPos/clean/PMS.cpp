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

	LevelLimit=6;
	paOrders=genPAorders(LevelLimit);
	
}


PMStruc::PMStruc(PyrMode p, int lvlimt)
{

	mymode=p;
	numOfData=0;

	LevelLimit= lvlimt;
	paOrders=genPAorders(LevelLimit);
	
	
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




int PMStruc::generatePymFromdata(vector<vector<double> > data,int adm)
{

	adimension=adm;
	return dataToPym(data);

}

vector<pair<int,int> > PMStruc:: dataToTwoInxs(vector<double> data)
{


	vector<pair<int,int> > inxs;
	inxs.clear();
	for (int i = 0; i <LevelLimit ; i++)
	{
		inxs.push_back( dataToTwoInxWthoutPos(i,data) );
	}

	
	
	return inxs;
}

pair<int,int> PMStruc::dataToTwoInxWthoutPos(int alvel,vector<double> data)
{

	int dimension=data.size()+2;
	int haldim=dimension/2;

	int fisI(0),secI(0);


	
	for (int j=2;j<haldim;j++)
	{
		int indtem;
		indtem=aAbs[alvel][j].first*data[j-2]+aAbs[alvel][j].second;
		if (indtem<0)
		{
			indtem=0;
		}
		if (indtem>(pow2[alvel]-1))
		{
			indtem=pow2[alvel]-1;
		}
		fisI+=twoExs[j]*indtem;
	}
	for (int j=haldim;j<dimension;j++)
	{
		int indtem;
		indtem=aAbs[alvel][j].first*data[j-2]+aAbs[alvel][j].second;
		if (indtem<0)
		{
			indtem=0;
		}
		if (indtem>(pow2[alvel]-1))
		{
			indtem=pow2[alvel]-1;
		}
		secI+=twoExs[j]*indtem;
	}


	return pair<int,int>(fisI,secI);

	

		

}



pair<int,int> PMStruc::	lvelToAPlvel(int lvel)
{
	int alvel;
	int plvel;

	alvel=lvel%LevelLimit;
	
	plvel=lvel/LevelLimit;
	return pair<int,int>(alvel,plvel);
}

pair<int,int> PMStruc::dataToTwoPosInx(int lvel,vector<double> data,vector<pair<int,int> > inp)
{
	pair<int,int> apls=lvelToAPlvel(lvel);
	int alvel=apls.first;
	int plvel=apls.second;

	

	int dimension=data.size();
	int haldim=dimension/2;

	int fisI(0),secI(0);


	for (int j = 0; j < 2; j++)
	{
		int indtem;
		indtem=aAbs[plvel][j].first*data[j]+aAbs[plvel][j].second;
		if (indtem<0)
		{
			indtem=0;
		}
		if (indtem>(twExMs-1))
		{
			indtem=twExMs-1;
		}
		fisI+=twoExs[j]*indtem;
	}
	return pair<int,int>(inp[alvel].first+fisI,inp[alvel].second);
}

pair<int,int> PMStruc::dataToTwoPosInx(int lvel,vector<double> data)
{
	pair<int,int> apls=lvelToAPlvel(lvel);
	int alvel=apls.first;
	int plvel=apls.second;

	

	int dimension=data.size();
	int haldim=dimension/2;

	int fisI(0),secI(0);


	for (int j = 0; j < 2; j++)
	{
		int indtem;
		indtem=aAbs[plvel][j].first*data[j]+aAbs[plvel][j].second;
		if (indtem<0)
		{
			indtem=0;
		}
		if (indtem>(twExMs-1))
		{
			indtem=twExMs-1;
		}
		fisI+=twoExs[j]*indtem;
	}
	for (int j=2;j<haldim;j++)
	{
		int indtem;
		indtem=aAbs[alvel][j].first*data[j]+aAbs[alvel][j].second;
		if (indtem<0)
		{
			indtem=0;
		}
		if (indtem>(twExMs-1))
		{
			indtem=twExMs-1;
		}
		fisI+=twoExs[j]*indtem;
	}
	for (int j=haldim;j<dimension;j++)
	{
		int indtem;
		indtem=aAbs[alvel][j].first*data[j]+aAbs[alvel][j].second;
		if (indtem<0)
		{
			indtem=0;
		}
		if (indtem>(twExMs-1))
		{
			indtem=twExMs-1;
		}
		secI+=twoExs[j]*indtem;
	}


	return pair<int,int>(fisI,secI);
}

/*
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
*/

void valueToInx(pair<double,double> minMax,pair<double,double>& aAndB,int levl)
{
	int totl=pow2[levl];
	minMax.first-=0.001;
	minMax.second+=0.001;
	aAndB.first=(double)(totl)/(minMax.second-minMax.first);
	aAndB.second=0-aAndB.first*minMax.first;
}



bool PMStruc::dataToPymLvl(vector<vector<double> > datas,int lvel)
{
	bool fineEngough=true;
	int dimension=datas[0].size();
	
	int haldim=datas[0].size()/2;


	for (int i=0;i<datas.size();i++)
	{
		pair<int,int> ss;
	
		ss= dataToTwoPosInx(lvel,datas[i]);
		int fisI(ss.first),secI(ss.second);

		if (pym[lvel].count(fisI)>0)
		{
			if (pym[lvel][fisI].count(secI)>0)
			{
				pym[lvel][fisI][secI]+=1;
				fineEngough=false;
			}
			else
			{
				pym[lvel][fisI].insert(pair<int,int>(secI,1));
			}
		}
		else
		{
			unordered_map<int,int> temlv;
			temlv.insert(pair<int,int>(secI,1));
			pym[lvel].insert(pair<int,unordered_map<int,int> >(fisI,temlv));
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
/*
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
*/
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
		
		twExMs=pow2[LevelLimit-1];
		twoExs.resize(dimension,0);
		for (int ti=0;ti<haldim;ti++)
		{
			twoExs[ti]=twoExs[ti+haldim]=pow(twExMs,ti);
		}
		for (int goodi = 0; goodi < LevelLimit; goodi++)
		{


			vector<pair<double,double> > aAb;
			aAb.resize(data[0].size(),pair<double,double>(0.0,0.0));

			for (int i=0;i<data[0].size();i++)
			{
				valueToInx(minmax[i],aAb[i],goodi);
			}

			aAbs.push_back(aAb);
			
	

		}

	
		weights.resize(LevelLimit*LevelLimit,0.0);
		unordered_map<int,unordered_map<int,int>> tum;
		pym.resize(LevelLimit*LevelLimit,tum);
		
		for (int i = 0; i < LevelLimit*LevelLimit; i++)
		{



				dataToPymLvl(data,i);
				pair<int,int> apls=lvelToAPlvel(i);

				weights[i]=sqrt((double)(pow2[apls.first]*pow2[apls.second]));


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

/*vector<int> PMStruc:: matchDToOneLv(vector<vector<double> > dataset,int levl,map<int,map<int,int> > pmlv,bool ExcluMode,set<int> elimitset )
{

	int haldim=dataset[0].size()/2;
	int dimension=dataset[0].size();

	int siz=dataset.size();
	vector<int> res;
	res.clear();
	for (int i=siz-1;i>=0;i--)
	{
		
		pair<int,int> ss;
	
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
					//	res+=1;
						pmlv[fisI][secI]-=1;
						if (elimitset.count(i)==0)
						{
							res.push_back(i);
						}
					//	dataset.erase(dataset.begin()+i);
					}
				}
				else 
					{
						if (pmlv[fisI][secI]>0)
						{
						//	res+=1;
							pmlv[fisI][secI]-=1;
							if (elimitset.count(i)==0)
							{
								res.push_back(i);
							}
						//	dataset.erase(dataset.begin()+i);
						}
				}
			}
		}
	}
	return res;
}
*/


int PMStruc:: matchDToOneLvSimple(vector<vector<pair<int,int> > > vecFSs,vector<vector<double>> poss,int levl, bool ExcluMode,vector<bool>& used )
{

	


	int siz=vecFSs.size();
	int res(0);

	auto pmlv=pym[levl];

	for (int i=siz-1;i>=0;i--)
	{
		if(!used[i])
		{
			pair<int,int> ss;
	
		ss=dataToTwoPosInx(levl,poss[i],vecFSs[i]);

		int fisI(ss.first),secI(ss.second);

		if (pmlv.count(fisI)>0)
		{
			if (pmlv[fisI].count(secI)>0)
			{
				if(ExcluMode)
				{
					if (pmlv[fisI][secI]>1)
					{
						
						if (used[i]==false)
						{
							res+=1;
							used[i]=true;
						}
					}
				}
				else 
					{
						if (pmlv[fisI][secI]>0)
						{
							if (used[i]==false)
							{
								res+=1;
								used[i]=true;
							}
						
						}
				}
			}
		}
		}
		
	}
	return res;


	
	return 0;
}

double PMStruc::MatchDttoPosPymSimple(vector<vector<pair<int,int> > > vecFSs,vector<vector<double>> poss,bool ExcluMode,vector<double> & mnumbers,bool inverse)
{
	int siz=vecFSs.size();

	
	mnumbers.resize(pym.size(),0.0);

	vector<bool> used;
	used.clear();
	used.resize(vecFSs.size(),false);


	for (int i = 0; i < paOrders.size(); i++)
	{
			for (int j = 0; j < paOrders[i].size(); j++)
			{
				int tinx=paOrders[i][j].first*LevelLimit+paOrders[i][j].second;
				mnumbers[tinx]=matchDToOneLvSimple(vecFSs,poss,tinx,ExcluMode,used);
			}
	}
	

	double reslt(0.0);
	for (int i=0;i<pym.size();i++)
	{
		if(!inverse)
			reslt+=mnumbers[i]*weights[i]*(adimension+2);
		else
			reslt+=mnumbers[i]*(1.0/weights[i])*(adimension+2);
	}

	reslt/=siz;
	if(inverse)
	{
		
		reslt=1.0/reslt;
	}

		
	return reslt;
	
	return 0.0;
}






double PMStruc::MatchDttoPosPymSimple(vector<vector<double> > dataset,bool ExcluMode,vector<double> & mnumbers,bool inverse)
{
	auto apPs=toTwo(dataset,2);

	vector<vector<pair<int,int> > > faSs;
	//faSs.resize(dataset.size(),vector<pair<int,int> > ());
	faSs.clear();
	for (int i = 0; i < apPs.second.size(); i++)
	{
		faSs.push_back( dataToTwoInxs (apPs.second[i]) );
	}

	

		
	return MatchDttoPosPymSimple(faSs,apPs.first,ExcluMode,mnumbers,inverse);
}
/*double PMStruc::MatchDttoPosPym(vector<vector<double> > dataset,bool ExcluMode,vector<double> & mnumbers,bool inverse)
{
	set<int> used;
	used.clear();
	mnumbers.resize(pym.size(),0.0);
	double reslt(0.0);
	for (int i = 0; i < paOrders.size(); i++)
	{
		vector<vector<int> > lrslt;
		lrslt.resize(paOrders[i].size(),vector<int>());


		for (int j = 0; j < paOrders[i].size(); j++)
		{
			int tinx=paOrders[i][j].first*LevelLimit+paOrders[i][j].second;
			lrslt[j]=matchDToOneLv(dataset,tinx,pym[tinx],ExcluMode,used);
		}
		
		vector<vector<int> > belongTo;
		belongTo.resize(dataset.size(),vector<int>());
		
		for (int j = 0; j < lrslt.size(); j++)
		{
			for (int k = 0; k < lrslt[j].size(); k++)
			{
				belongTo[lrslt[j][k]].push_back(j);
			}
		}

		for (int j = 0; j < belongTo.size(); j++)
		{
			if (belongTo[j].size()>0)
			{
				for (int k = 0; k < belongTo[j].size(); k++)
				{
					int mj=belongTo[j][k];
					int tinx=paOrders[i][mj].first*LevelLimit+paOrders[i][mj].second;
					mnumbers[tinx]+=1.0/(double)belongTo[j].size();
				}
				used.insert(j);
			}
		}
	}
	int dim=dataset[0].size();
	for (int i = 0; i < mnumbers.size(); i++)
	{
		if(!inverse)
			reslt+=mnumbers[i]*weights[i]*dim;
		else
			reslt+=mnumbers[i]*(1.0/weights[i])*dim;
	}
	int siz=dataset.size();
	reslt/=siz;
	if(inverse)
	{
		
		reslt=1.0/reslt;
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

void setDoubleVec(vector<double>& to,vector<double> frm)
{
	to.clear();
	to.insert(to.end(),frm.begin(),frm.end());
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

PMSEnsemble::PMSEnsemble()
{
	pyms.clear();
}

int PMSEnsemble::generateAaBsFromdata(vector<vector<double> > data,int LevelLimit)
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
	
		twExMs=pow2[LevelLimit-1];
		twoExs.resize(dimension,0);
		for (int ti=0;ti<haldim;ti++)
		{
			twoExs[ti]=twoExs[ti+haldim]=pow(twExMs,ti);
		}

		weights.resize(LevelLimit*LevelLimit,0.0);
		unordered_map<int,unordered_map<int,int>> tum;
		pym.resize(LevelLimit*LevelLimit,tum);
		
		for (int i = 0; i < LevelLimit*LevelLimit; i++)
		{
				int ai(i/LevelLimit);int aj(i%LevelLimit);
				weights[i]=sqrt((double)(pow2[ai]*pow2[aj]));
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
	
		pair<int,int> ss;
	
		ss= dataToTwoPosInx(lvel,data);
		int fisI(ss.first),secI(ss.second);
		

		
		if (pym[levi].count(fisI)>0)
		{
			if (pym[levi][fisI].count(secI)>0)
			{
				if(AddOrMinus)
					pym[levi][fisI][secI]+=1;
				else
					pym[levi][fisI][secI]-=1;
				
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
			unordered_map<int,int> temlv;
			if(AddOrMinus)
				temlv.insert(pair<int,int>(secI,1));
			else
				temlv.insert(pair<int,int>(secI,-1));
			pym[levi].insert(pair<int,unordered_map<int,int> >(fisI,temlv));
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
				PMStruc tp(PMStruc::postitionSpecific);
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
					vector<double> temiv;
					tr[i]=pyms[i].givePyramidMatchScore(d,false,temiv)/sqrt((double)pyms[i].getNumofData());
					
				}
				auto sd=maxminValAndInx(tr,true);
				if(sd.second>threshold)
					pyms[sd.first].AddSeverlData(d,true);
				else
				{
					
					PMStruc tp(PMStruc::postitionSpecific);
					tp.initPymWithABs(aAbs,dimension);
					tp.AddSeverlData(d,true);
					pyms.push_back(tp);
				}


			}

	

	return 0;
}

double PMStruc::givePyramidMatchScore(vector<vector<double> > dataset,bool ExcluMode,vector<double>& scoreAllLevel)
{

	return MatchDttoPosPymSimple(dataset,ExcluMode,scoreAllLevel,false);

}

double PMStruc::givePyramidMatchScore(vector<vector<pair<int,int> > > vecFSs,vector<vector<double>> poss,bool ExcluMode,vector<double>& scoreAllLevel)
{
	assert(vecFSs.size()==poss.size());
//double MatchDttoPosPymSimple(vector<pair<int,int>> vecFSs,vector<vector<double>> poss,bool ExcluMode,vector<double> & mnumbers,bool inverse);
	return  MatchDttoPosPymSimple(vecFSs,poss,ExcluMode,scoreAllLevel,true);
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
			rslts.push_back(s.givePyramidMatchScore(dataset,false,vector<double>()));

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