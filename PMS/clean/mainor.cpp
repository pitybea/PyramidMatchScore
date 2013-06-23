/*
 * main.cpp
 *
 *  Created on: Jul 13, 2012
 *      Author: pitybea
 */

#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <math.h>
#include <stdio.h>


#include "PMS.h"
#include "codeSpe.h"

using namespace std;
using namespace cv;




void getMatric(vector<vector<int> >&vls, string fname)
{
	ifstream iof;
	iof.open(fname.c_str());
	int wid,hei;
	iof>>wid>>hei;
	vls.clear();
	vls.resize(hei,vector<int>(wid,0));
	for (int i=0;i<hei;i++)
	{
		for (int j=0;j<wid;j++)
		{
			iof>>vls[i][j];
		}

	}

}
void ptsTov(vector<Point> p,vector<vector<double> >& d)
{
	d.resize(p.size(),vector<double>(2,0.0));
	for (int i=0;i<p.size();i++)
	{
		d[i][0]=p[i].x;
		d[i][1]=p[i].y;
	}

}

int main()
{



	vector<vector<int> > tmatch;
	vector<vector<int> > templ;
	getMatric(templ,"m1.txt");
	getMatric(tmatch,"m2.txt");

	vector<vector<double> > dats;

	for (int i=0;i<templ.size();i++)
	{
		for (int j=0;j<templ[0].size();j++)
		{
			if (templ[i][j]>0)
			{
				vector<double> tddd;
				tddd.push_back(j*1.00);
				tddd.push_back(i*1.00);
				dats.push_back(tddd);
			}

		}

	}


	vector<Point> egpts;
	egpts.clear();
	for (int i=0;i<tmatch.size();i++)
	{
		for (int j=0;j<tmatch[0].size();j++)
		{
			if (tmatch[i][j]>0)
			{
				Point p;
				p.y=i;
				p.x=j;
				egpts.push_back(p);
			}

		}

	}





//	int totalLvls(4);

	PMStruc pedpm(4);


	pedpm.dataToPym(dats);

	vector<Rect> recs;
	recs.clear();
	vector<double> scores;
	scores.clear();





	int curx(0),cury(0);
	double scale(1.0);
	double scaleInve=1/scale;

	int rwidth(100),rhei(40);
	int stepx(2),stepy(2);

	int imgw=640;
	int imgh=480;


	rwidth*=scale;
	rhei*=scale;

	cout<<cury<<" "<<rhei<<endl;

	while ((cury+rhei)<tmatch.size()+40)
	{
		curx=0;
		cout<<cury<<endl;
		while((curx+rwidth)<tmatch[0].size()+100)
		{


			vector<Point> theMpss;
			theMpss.clear();
			for (int i=0;i<egpts.size();i++)
			{
				int cex=curx+rwidth/2;
				int cey=cury+rhei/2;
				if (Rect(curx,cury,rwidth,rhei).contains(egpts[i]))
				{
					theMpss.push_back(Point((cex-egpts[i].x)*scaleInve,(cey-egpts[i].y)*scaleInve)) ;

				}
			}
			vector<vector<double> > theMsetP;
			ptsTov(theMpss,theMsetP);

			if(theMsetP.size()>0)
			{


				double temscore=pedpm.MatchDttoPym(theMsetP);
				{
					scores.push_back(temscore);
					recs.push_back(Rect(curx,cury,rwidth,rhei));
				}


			}

			curx+=stepx;
		}
		cury+=stepy;
	}





	int totalFod(0);
	vector<Rect> fRecs;
	vector<double> fScors;
	fRecs.clear();
	fScors.clear();
	while((totalFod<10) && (scores.size()>0))
	{
		int maXi=MAxInd(scores);
		fScors.push_back(scores[maXi]);
		fRecs.push_back(recs[maXi]);


		for (int i=egpts.size()-1;i>=0;i--)
		{

			if (recs[maXi].contains(egpts[i]))
			{
				egpts.erase(egpts.begin()+i);
			}

		}
		recs.erase(recs.begin()+maXi);

		totalFod+=1;
		scores.clear();
		vector<Rect> temRcs;
		temRcs.clear();
		temRcs.swap(recs);
		recs.clear();

		for (int ri=0;ri<temRcs.size();ri++)
		{

			vector<Point> theMpss;
			theMpss.clear();
			for (int i=0;i<egpts.size();i++)
			{
				if (temRcs[ri].contains(egpts[i]))
				{
					theMpss.push_back(Point(temRcs[ri].x+temRcs[ri].width/2-egpts[i].x
						,temRcs[ri].y+temRcs[ri].height/2-egpts[i].y));
				}

			}
			vector<vector<double> > theMsetP;
			ptsTov(theMpss,theMsetP);
			if(theMsetP.size()>0)
			{

				double temscore=pedpm.MatchDttoPym(theMsetP);//-backpm.MatchDttoPym(theMset);
//				if (temscore>100)
				{
					scores.push_back(temscore);
					recs.push_back(temRcs[ri]);
				}


			}
		}



	}


		scores.swap(fScors);
		recs.swap(fRecs);

		vector<int> inde_;
		inde_.resize(scores.size(),0);
		for (int i=0;i<scores.size();i++)
		{
			inde_[i]=i;
		}
		prshl(scores,scores.size(),inde_);

		FILE* ftt;
		string tfn="_pms_s8_LaRec";
		ftt=fopen(tfn.c_str(),"w");
		fprintf(ftt,"%d\n",scores.size());
		for (int i=0;i<40 && i<scores.size();i++)
		{
			int j=inde_[inde_.size()-1-i];
			//rectangle(iptimg,recs[j],Scalar(255,255,0),1,8,0);
			fprintf(ftt,"%d %d %d %d %lf\n",recs[j].x,recs[j].y,recs[j].width,recs[j].height,scores[j]);
		}
		string savInl;
//		savInl=tD+"_s10_rec_po.jpg";
	//	imwrite(savInl,iptimg);
		fclose(ftt);
		cout<<endl;




	return 0;

}


