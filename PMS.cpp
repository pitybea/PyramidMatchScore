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
	name="";
}
PMStruc::PMStruc(int i,PyrMode p,string s)
{
	totalLvls=i;
	mymode=p;
	name=s;
	weights.resize(i,0.0);
	for(int i=0;i<weights.size();i++)
	{
		weights[i]=pow2[i];
	}
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
	fprintf(fp,"%s %d %d\n%d\n",name,mymode,totalLvls,pym.size());
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
	string wFname=name+"_wgt.txt";
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
	fscanf(fp,"%s %d %d\n%d\n",&tems,&mymode,&totalLvls,&temsz);
	string s(tems);
	name=s;
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
	string wFname=name+"_wgt.txt";
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
