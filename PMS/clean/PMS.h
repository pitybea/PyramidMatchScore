#include <vector>
#include <string>
#include <map>
#include <math.h>
#include <assert.h>
using namespace std;



#include <set>

struct Point
{
	int x;
	int y;
};

template<class T>
vector<vector<T> > selectVecButLstTwo(vector<vector<T> > inp,int dim)
{
	
	vector<vector<T> > rslt;
	rslt.clear();

	if(inp.size()>0)
	{
		auto myd=inp[0].size();
		for (auto ss:inp)
		{
			 ss.erase(ss.begin()+dim, ss.begin()+myd-2 );
			rslt.push_back(ss );
		}
	}
	return rslt;
}

template<class T>
vector<vector<T> > selectVec(vector<vector<T> > inp,int dim)
{
	
	vector<vector<T> > rslt;
	rslt.clear();

	if(inp.size()>0)
	{
		auto myd=inp[0].size();
		for (auto ss:inp)
		{
			 ss.erase(ss.begin()+dim, ss.begin()+myd );
			rslt.push_back(ss );
		}
	}
	return rslt;
}


template<class T>
static void prshl(vector<T> p,int n,vector<int>& index)
{
	int k,j,i;
	T t;
	int ii;
	k=n/2;
	while(k>0)
	{
		for(j=k;j<=n-1;j++)
		{
			t=p[j];  ii=index[j];  i=j-k;
			while((i>=0)&&(p[i]>t))
			{
				p[i+k]=p[i];  index[i+k]=index[i];  i=i-k;
			}
			p[i+k]=t;  index[i+k]=ii;
		}
		k=k/2;
	}
};



template <class T,class U>
int bendPsnToApea(vector<vector<T> > &features,vector<U> pts)
{
	if (features.size()==pts.size())
	{
		for (int i=0;i<pts.size();i++)
		{
			features[i].push_back(pts[i].x);
			features[i].push_back(pts[i].y);
		}
	}
	return 0;
}

template<class T>
int NormalVec(vector<vector<T> >& dataset)
{
	if (dataset.size()>0)
	{
		if (dataset[0].size()>0)
		{
			for (int i=0;i<dataset.size();i++)
			{
				double temtot(0.0);
				for (int j=0;j<temtot;j++)
				{
					temtot+=dataset[i][j];
				}
				temtot+=0.0001;
				for (int j=0;j<temtot;j++)
				{
					dataset[i][j]/=temtot;
				}
			}
		}
	}
	return 0;
}



template<class T>
vector<vector<T> > TransitMtx(vector<vector<T> > oData,vector<vector<T> >trasM)
{
	vector<vector<T> > tData;
	tData.resize(oData.size(),vector<T>(trasM[0].size(),0.0));
	
	for (int i=0;i<oData.size();i++)
	{

		for (int j=0;j<trasM[0].size();j++)
		{
			for (int k=0;k<trasM.size();k++)
			{
				tData[i][j]+=oData[i][k]*trasM[k][j];
			}
		}
	}
	return tData;
}


template<class T>
pair<int,T> maxminValAndInx(vector<T> inp,bool maxormin)
{
	pair<int,T> rslt;
	if(inp.size()==0)
	{
		rslt.first=-1;
		return rslt;
	}
	else if(inp.size()==1)
	{
		rslt.first=0;
		rslt.second=inp[0];
		return rslt;

	}
	else
	{
		rslt.first=0;
		rslt.second=inp[0];
		for (int i = 1; i < inp.size(); i++)
		{
			if(maxormin)
			{
				if (inp[i]>rslt.second)
				{
					rslt.first=i;
					rslt.second=inp[i];
				}
			}
			else
			{
				if (inp[i]<rslt.second)
				{
					rslt.first=i;
					rslt.second=inp[i];
				}
			}
		}

		return rslt;
	}

	

}

template<class T>
int ZeroMnVec(vector<vector<T> >& dataset)
{
	if (dataset.size()>0)
	{
		if (dataset[0].size()>0)
		{
			for (int i=0;i<dataset.size();i++)
			{
				double temtol(0.0);
				for (int j=0;j<dataset[i].size();j++)
				{
					temtol+=dataset[i][j];
				}
				temtol/=dataset[0].size();
				for (int j=0;j<dataset[i].size();j++)
				{
					dataset[i][j]-=temtol;
				}
			}
		}
	}
	return 0;
}

template <class T>
int MAxInd(vector<T> data)
{
	T a;
	a=data[0];
	int ind=0;
	for (int i=0;i<data.size();i++)
	{
		if (data[i]>a)
		{
			a=data[i];
			ind=i;
		}
	}
	return ind;
};


struct featype
{
	featype()
	{

	};
	Point pos;
	double feature[128];
	static void initOne(FILE* fp,featype &t)
	{
		fscanf(fp,"%d %d\n",&t.pos.x,&t.pos.y);
		for (int i=0;i<128;i++)
		{
			fscanf(fp,"%lf ",&t.feature[i]);
		}
	};
	static void printOne(FILE* fp, featype t)
	{
		fprintf(fp,"%d %d\n",t.pos.x,t.pos.y);
		for (int i=0;i<128;i++)
		{
			fprintf(fp,"%lf ",t.feature[i]);
		}
	};
	vector<double> toVdouble()
	{
		vector<double> res;
		res.resize(128,0.0);
		for (int i=0;i<128;i++)
		{
			res[i]=feature[i];
		}
	//	res.push_back(pos.x);
	//	res.push_back(pos.y);

		return res;
	};

};

static vector<vector<double> > keepFirSevDims(vector<vector<double> > inp,int n)
{
	vector<vector<double> > results;
	results.resize(inp.size(),vector<double>(0,0.0));
	for (int i = 0; i < inp.size(); i++)
	{
		results[i].insert(results[i].end(),inp[i].begin(),inp[i].begin()+n);
	}
	return results;
}

static vector<vector<double> > addPositionsToData(vector<vector<double> > data,vector<featype> feas)
{
	vector<vector<double> > result(data);
	//	result=;
	assert(data.size()==feas.size());
	for (int i=0;i<result.size();i++)
	{
		result[i].push_back(feas[i].pos.x);
		result[i].push_back(feas[i].pos.y);
	}
	return result;
}

static int posorder[36][2]=
	{
	5,5,

	4,5,
	5,4,

	5,3,
	4,4,
	3,5,

	2,5,
	3,4,
	4,3,
	5,2,

	5,1,
	4,2,
	3,3,
	2,4,
	1,5,

	0,5,
	1,4,
	2,3,
	3,2,
	4,1,
	5,0,

	4,0,
	3,1,
	2,2,
	1,3,
	0,4,

	0,3,
	1,2,
	2,1,
	3,0,

	2,0,
	1,1,
	0,2,

	0,1,
	1,0,

	0,0
	};

static vector<vector<pair<int,int> > > genPAorders(int n)
{
	int m=2*(n-1);
	vector<vector<pair<int,int> > > rslt;

	for (int i = m; i >=0; i--)
	{
		vector<pair<int,int> > lvl;
		lvl.clear();

		int s=((i>=(n-1))?(n-1):i);
		int t=i-s;

		while ( (t>=0) && (s<=(n-1)) && (t<=(n-1)) && (s>=0) )
		{
			lvl.push_back(pair<int,int>(s,t));
			s-=1;
			t=i-s;
		}


		rslt.push_back(lvl);

	}
	return rslt;
};

class PMStruc
{
public:
	enum PyrMode{normal,/*average,*/inverse,postitionSpecific};
	PMStruc(PyrMode p);
	PMStruc(PyrMode p,int);
	PMStruc(){};

	PyrMode mymode;
	//string name;

	vector<double> weights;

	 
	int generatePymFromdata(vector<vector<double> > data);

	//int generatePosPymsFromdata(vector<vector<double> > data);


	double givePyramidMatchScore(vector<vector<double> > dataset,bool ExcluMode,vector<double>& scoreAllLevel);
	
	void outToAFile(string filename);
	void loadFromAFile(string filename);
	

	
	int initPymWithABs(vector<vector<pair<double,double> > > abS,int dimension);
	int AddoneData(vector<double> data,bool AddOrMinus);
	int AddSeverlData(vector<vector<double> > data,bool AddorMinus);
	int getNumofData(){return numOfData;};
	int GeneratePosWeightWithParameter(double ratio);


private:
	int dataToPym(vector<vector<double> > data);
//	int dataToPymAver(vector<vector<double> > data);
//	int dataToPosPyms(vector<vector<double> > data);

	//void valueToInx(pair<double,double> minMax,pair<double,double>& aAndB,int levl);

	bool dataToPymLvl(vector<vector<double> > datas,int lvel,map<int,map<int,int> >& pymlvl);
//	bool dataToPymLvl(vector<vector<double> > datas,int lvel,map<int,map<int,int> >& pymlvl,vector<vector<double> > aintvl);
//	bool dataToPosPymLvl(vector<vector<double> > datas,int alvel,int plvel,map<int,map<int,int> >& pymlvl);


	double MatchDttoPym(vector<vector<double> > dataset,bool ExcluMode,vector<double> & scoreAllLevel,bool inverse);
//	double MatchDttoPymAv(vector<vector<double> > dataset,bool ExcluMode,vector<int> & scoreAllLevel);
//	double MatchPosDttoPym(vector<vector<double> > dataset,bool ExcluMode,vector<int> & mnumbers,bool inverse,int order[LevelLimit*LevelLimit][2]);
	double MatchDttoPosPym(vector<vector<double> > dataset,bool ExcluMode,vector<double> & mnumbers,bool inverse);

	int matchDToOneLv(vector<vector<double> > dataset,int levl,map<int,map<int,int> > pmlv, bool ExcluMode );

	vector<int> matchDToOneLv(vector<vector<double> > dataset,int levl,map<int,map<int,int> > pmlv,bool ExcluMode,set<int> elimitset );
//	int matchDToOneLv(vector<vector<double> > & dataset,int levl,map<int,map<int,int> > pmlv,vector<vector<double> > invs, bool ExcluMode );
//	int matchDToOnePosLv(vector<vector<double> >& dataset,int alevl,int plevl,map<int,map<int,int> > pmlv, bool ExcluMode );
	

	pair<int,int> dataToTwoIndx(int,vector<double>);
	pair<int,int> dataToTwoPosInx(int alvel,int plvel,vector<double> data);

	vector<map<int,map<int,int> > > pym;

	//vector<vector<map<int,map<int,int> > > > pospyms;

	vector<vector<pair<double,double> > > aAbs;
//	vector<vector<vector<double> > > intvDecs;

	//static int posSpcOrder[36][2](posorder);
	
	//int haldim;
	vector<int> twExMs;//=pow2[i];
	vector<vector<int>> twoExs;

	int numOfData;
	int LevelLimit;// 6
	vector<vector<pair<int,int> > > paOrders;
};

class PMSEnsemble
{
public:
	PMSEnsemble();
	vector<vector<double>> weights;
	int generateAaBsFromdata(vector<vector<double> > data,int);

	int generateStructureFromData(vector<vector<vector<double> > > data);
	double givePyramidMatchScore(vector<vector<double> > dataset);
	
	vector<vector<pair<double,double> > > aAbs;
	double threshold;


//private:

	vector<PMStruc> pyms;
	int dimension;
	
};

/*
class PMStrainer
{
public:
	enum TrainStrategy
	{
		allAtOnce,
		oneAtOnce

	};
	PMStrainer(PMStruc pm,TrainStrategy s,double d,double l,string ts);
	void Train(vector<vector<int> > posnums,vector<vector<int> > negnums);
	void OutTofile();
	string name;



private:
	PMStruc pmStrct;
	TrainStrategy myStrgy;
	double stepWidth;
	double largeBand;
	vector<double> bestWeights;
	double bestDivScore;
	int bestNumber;

	void trainOneAtOnce(vector<vector<int> > posnums,vector<vector<int> > negnums);
	void trainAllAtOnce(vector<vector<int> > posnums,vector<vector<int> > negnums);
	int DecideBestScore(vector<double> temWeight,vector<vector<int> > posnums,vector<vector<int> > negnums,double& divValue);
};
*/

class trainVecs
{
public:
	enum strategies{allAtOnce,oneAtonce};

	trainVecs(strategies strg,vector<vector<double> > lhbds,vector<double> stps,vector<vector<double> > pexs,vector<vector<double> > nexs);

	void train();
	void printWeights()
	{
		for(auto ss:bstWeights) 
		printf("%lf ",ss);
	};
private:
	strategies mystr;
	vector<vector<double> > posvecs;
	vector<vector<double> > negvecs;
	vector<double> bstWeights;
	vector<double> stps;
	vector<vector<double> > lowHighBands;

	void trainOneAtOnce();
	void trainAllAtOnce();
};


static vector<vector<double > > prepareData(vector<featype> allfeas, vector<vector<double> > TrmTx, bool addPos)
{
	vector<vector<double > >  dats;
	dats.resize(allfeas.size(),vector<double>(0,0.0));
	for (int i=0;i<dats.size();i++)
	{
		dats[i]=allfeas[i].toVdouble();
	}
	
	ZeroMnVec(dats);
	dats=TransitMtx(dats,TrmTx);
	if(addPos)
		dats=addPositionsToData(dats,allfeas);
	return dats;
}

