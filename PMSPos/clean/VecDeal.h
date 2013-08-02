#include <vector>
#include <string>
#include <map>
#include <math.h>
#include <assert.h>
#include <unordered_map>

using namespace std;



#include <set>

struct Point
{
	int x;
	int y;
};

//this function will put position information at the begining of a vector (which is the last two dimensions)
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
			vector<T> rvec;
			rvec.clear();
			rvec.insert(rvec.end(),ss.begin()+myd-2,ss.end());
			ss.erase(ss.begin()+dim, ss.end());

			rvec.insert(rvec.end(),ss.begin(),ss.end());
			

			rslt.push_back(rvec );
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

template<class T>
static pair<vector<vector<T> >,vector<vector<T> > > toTwo(vector<vector<T> > inp, int dim)
{
	pair<vector<vector<T> >,vector<vector<T> > > rslt;


	for(auto ss: inp)
	{
		rslt.first.push_back(vector<T>(ss.begin(),ss.begin()+dim));
		rslt.second.push_back(vector<T>(ss.begin()+dim,ss.end()) );
	}
	return rslt;
}