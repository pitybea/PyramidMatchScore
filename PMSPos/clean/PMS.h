#include "VecDeal.h"





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
	enum PyrMode{postitionSpecific,positionsimple};
	PMStruc(PyrMode p);
	PMStruc(PyrMode p,int);
	PMStruc(PyrMode p,int,int);
	PMStruc(){};

	PyrMode mymode;
	//string name;

	vector<double> weights;
	 
	int generatePymFromdata(vector<vector<double> > data);
	
	double givePyramidMatchScore(vector<vector<double> > dataset,bool ExcluMode,vector<double>& scoreAllLevel);

	double givePyramidMatchScore(vector<vector<pair<int,int> > > vecFSs,vector<vector<double>> poss,bool ExcluMode,vector<double>& scoreAllLevel);


	int initPymWithABs(vector<vector<pair<double,double> > > abS,int dimension);
	int AddoneData(vector<double> data,bool AddOrMinus);
	int AddSeverlData(vector<vector<double> > data,bool AddorMinus);
	int getNumofData(){return numOfData;};
	

	
	int GeneratePosWeightWithParameter(double ratio);

	
	vector<pair<int,int> > dataToTwoInxs(vector<double> data);



private:
	int dataToPym(vector<vector<double> > data);

	bool dataToPymLvl(vector<vector<double> > datas,int lvel);



//	double MatchDttoPosPym(vector<vector<double> > dataset,bool ExcluMode,vector<double> & mnumbers,bool inverse);

	double MatchDttoPosPymSimple(vector<vector<double> > dataset,bool ExcluMode,vector<double> & mnumbers,bool inverse);

	double MatchDttoPosPymSimple(vector<vector<pair<int,int> > > vecFSs,vector<vector<double>> poss,bool ExcluMode,vector<double> & mnumbers,bool inverse);

	int matchDToOneLvSimple(vector<vector<double> > dataset,int levl, bool ExcluMode,vector<bool>&  used);

	int matchDToOneLvSimple(vector<vector<pair<int,int> > > vecFSs,vector<vector<double>> poss,int levl, bool ExcluMode,vector<bool>& used );

//	vector<int> matchDToOneLv(vector<vector<double> > dataset,int levl,map<int,map<int,int> > pmlv,bool ExcluMode,set<int> elimitset );

	pair<int,int>  dataToTwoInx(int alvel,vector<double> data);
	pair<int,int> dataToTwoPosInx(int lvel,vector<double> data);

	void addpostoFS(pair<int,int>& fs,int lvel,vector<double> data);
	
	vector<unordered_map<int,unordered_map<int,int> > > pym;



	vector<vector<pair<double,double> > > aAbs;

	


	pair<int,int> lvelToAPlvel(int);


	int twExMs;//=pow2[i];
	vector<int> twoExs;

	int adimension;
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

