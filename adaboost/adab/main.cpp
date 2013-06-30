#include "adaboost.h"
#include <Windows.h>
#include <direct.h>
#define dim_S_n 6
#define num_Tr_ePs 1050
#define num_Tes_ePs 889
#define num_Wk_cLsfs 5

adaboost ada;

void init()
{
	ada.dimension=dim_S_n;
	ada.num_smps=num_Tr_ePs;
	ada.num_P_smps=0;
	ada.num_N_smps=0;
	ada.num_wk_cl=num_Wk_cLsfs;

	vector<double> tmv_;
	tmv_.resize(dim_S_n,0.0);

	ada.train_label.clear();
	ada.train_smps.clear();

	ada.train_label.resize(num_Tr_ePs,0);
	ada.train_smps.resize(num_Tr_ePs,tmv_);

	FILE* fp;
	fopen_s(&fp,"allfeas.txt","r");
	for (int i=0;i<num_Tr_ePs;i++)
	{
		for (int j=0;j<dim_S_n;j++)
		{
			fscanf_s(fp,"%lf ",&ada.train_smps[i][j]);
		}
		fscanf_s(fp,"\n");
	}
	fclose(fp);

	fopen_s(&fp,"lbs.txt","r");

	for (int i=0;i<num_Tr_ePs;i++)
	{
		fscanf_s(fp,"%d\n",&ada.train_label[i]);
		if (ada.train_label[i]==1)
		{
			ada.num_P_smps+=1;
		}
		else
			ada.num_N_smps+=1;
	}
	fclose(fp);
}

int main(int argc,char* argv[])
{
	_chdir("E:\\CarData\\TrainImages");
	init();
	ada.training();
	FILE* fp;
//	fopen_s(&fp,"adabMach.txt","r");
//	ada.loadFromfile(fp);
//	fclose(fp);

	vector<double> tssm;
	tssm.resize(dim_S_n,0.0);

	
	fopen_s(&fp,"test.txt","r");
	for (int i=0;i<num_Tes_ePs;i++)
	{
		for (int j=0;j<dim_S_n;j++)
		{
			fscanf_s(fp,"%lf ",&tssm[j]);
		}
		fscanf_s(fp,"\n");
		printf("%d\n",ada.classfy(tssm));
	}
	fclose(fp);

	fopen_s(&fp,"adabMach.txt","w");
	ada.outputTofile(fp);
	fclose(fp);
	return 0;
}