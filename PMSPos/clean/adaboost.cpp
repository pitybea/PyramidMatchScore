#include "adaboost.h"

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


void adaboost::intitWts()
{
	if (train_smps.size()>0)
	{
		smp_weight.resize(train_smps.size(),0.0);
		double P_weight;
		double N_weight;
		P_weight=1/(2.0*num_P_smps);
		N_weight=1/(26.0*num_N_smps);
		for (unsigned int i=0;i<train_label.size();i++)
		{
			if (train_label[i]==1)
			{
				smp_weight[i]=P_weight;
			}
			else
				smp_weight[i]=N_weight;
		}
	}
}
void adaboost::ordering()
{
	orders.clear();
	if (train_smps.size()>0)
		if (dimension>0)
		{
			for (int i=0;i<dimension;i++)
			{
				vector<double> temD;
				vector<int> temI;
				temD.clear();
				temI.clear();
				temD.resize(train_smps.size(),0.0);
				temI.resize(train_smps.size(),0);
				for (unsigned int j=0;j<train_smps.size();j++)
				{
					temD[j]=train_smps[j][i];
					temI[j]=j;
				}
				prshl(temD,temD.size(),temI);
				orders.push_back(temI);

			}
		}
}

//direction:  ( >thresh => postive ) => ( dirction=1 )

double adaboost::train_wk_cl(int dim_)
{
	dircFisher temdF;
	temdF.dimensn=dim_;

	vector<double> sum_Perrors;
	vector<double> sum_Nerrors;
	sum_Nerrors.clear();
	sum_Nerrors.clear();
	sum_Perrors.resize(train_smps.size());
	sum_Nerrors.resize(train_smps.size());

	double totalPerrors(0.0);
	double totalNerrors(0.0);

	for (int i=0;i<train_smps.size();i++)
	{
		if(train_label[orders[dim_][i]]>0)
			totalPerrors+=smp_weight[orders[dim_][i]];
		else
			totalNerrors+=smp_weight[orders[dim_][i]];

		sum_Perrors[i]=totalPerrors;
		sum_Nerrors[i]=totalNerrors;
	}
	
	double min_error=-1.0;

	for (unsigned int i=0;i<train_smps.size()-1;i++)
	{
		double tem_thres=train_smps[orders[dim_][i]][dim_]+train_smps[orders[dim_][i+1]][dim_];
		double total_error;
		total_error=sum_Perrors[i]+(totalNerrors-sum_Nerrors[i]);
		if (min_error<0)
		{
			min_error=total_error;
			temdF.thresh=tem_thres;
			temdF.direction=1;

		}
		else if (total_error<min_error)
		{
			min_error=total_error;
			temdF.thresh=tem_thres;
			temdF.direction=1;

		}

		total_error=sum_Nerrors[i]+(totalPerrors-sum_Perrors[i]);;
		if (min_error<0)
		{
			min_error=total_error;
			temdF.thresh=tem_thres;
			temdF.direction=-1;

		}
		else if (total_error<min_error)
		{
			min_error=total_error;
			temdF.thresh=tem_thres;
			temdF.direction=-1;

		}
	}

	candida_str_cl.push_back(temdF);
	return min_error;

}
void adaboost::training()
{
	intitWts();
	ordering();

	str_cl.clear();
	//alphas.clear();

	vector<int> can_dim;
	can_dim.clear();

	for (int i=0;i<dimension;i++)
	{
		if (train_smps[orders[i][train_smps.size()-1]][i]>0.000001)
		{
			can_dim.push_back(i);
		}
	}
	
	double alp;
	while(str_cl.size()<num_wk_cl)
	{
		
		wk_cl_errors.clear();
		wk_cl_errors.resize(can_dim.size(),0.0);
		candida_str_cl.clear();

		for (unsigned int i=0;i<can_dim.size();i++)
		{
			wk_cl_errors[i]=train_wk_cl(can_dim[i]);
		}


		vector<int> in_;

		in_.clear();
		for (unsigned int i=0;i<wk_cl_errors.size();i++)
		{
			in_.push_back(i);
		}

		prshl(wk_cl_errors,wk_cl_errors.size(),in_);

		for (unsigned int i=0;i<can_dim.size();i++)
		{
			if (can_dim[i]==candida_str_cl[in_[0]].dimensn)
			{
//				can_dim.erase(can_dim.begin()+i);
				break;
			}
		}

		if (wk_cl_errors[in_[0]]<0.000001)
		{
			break;
		}
		else if (wk_cl_errors[in_[0]]>0.5)
		{
			num_wk_cl-=1;
			
			continue;
		}

		

		alp=0.5*log((1-wk_cl_errors[in_[0]])/wk_cl_errors[in_[0]]);
		candida_str_cl[in_[0]].alpha=alp;
		str_cl.push_back(candida_str_cl[in_[0]]);
	//	alphas.push_back(alp);

		double Zttem;
		Zttem=0.0;
		vector<double> zttts;
		zttts.clear();
		zttts.resize(train_smps.size(),0.0);

		for (unsigned int i=0;i<train_smps.size();i++)
		{
			zttts[i]=smp_weight[i]*exp(0.0-(str_cl[str_cl.size()-1].alpha*train_label[i]*(double)singleClassfy(str_cl[str_cl.size()-1],train_smps[i])));
			Zttem+=zttts[i];
		}
		for (unsigned int i=0;i<smp_weight.size();i++)
		{
			smp_weight[i]=zttts[i]/Zttem;
		}
		

	
	}
	
	


}

double adaboost::classfy(vector<double> vd)
{
	double totl(0.0);
	if (str_cl.size()>0)
		for (unsigned int i=0;i<str_cl.size();i++)
		{
			totl+=singleClassfy(str_cl[i],vd)*str_cl[i].alpha;
		}
	return totl;
}

int adaboost::singleClassfy(dircFisher dF,vector<double> smp)
{
	if (dF.direction==1)
	{
		if (smp[dF.dimensn]>dF.thresh)
		{
			return 1;
		}
		else
		{
			return -1;
		}
	}
	else
	{
		if (smp[dF.dimensn]>dF.thresh)
		{
			return -1;
		}
		else
		{
			return 1;
		}

	}
}
/*
double thresh;
int direction;
int dimensn;
double alpha;
*/

void adaboost::outputTofile(FILE* fp)
{
	fprintf_s(fp,"%d\n",str_cl.size());
	for (unsigned int i=0;i<str_cl.size();i++)
	{
		fprintf_s(fp,"%lf\t%d\t%d\t%lf\n",str_cl[i].thresh,str_cl[i].direction,str_cl[i].dimensn,str_cl[i].alpha);
	}
}

void adaboost::loadFromfile(FILE* fp)
{
	int temi;
	fscanf_s(fp,"%d\n",&temi);
	dircFisher temdF;
	str_cl.resize(temi,temdF);
	for (int i=0;i<temi;i++)
	{
		fscanf_s(fp,"%lf\t%d\t%d\t%lf\n",&str_cl[i].thresh,&str_cl[i].direction,&str_cl[i].dimensn,&str_cl[i].alpha);
	}
}