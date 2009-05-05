#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <assert.h>
#include "../include/3nodetri.h"
#include "../include/wTimer.h"

#define ELEMENT	node3tri//element type

using namespace std;
const int n=50;//
const double eps=1e-14;
const char* Resultsdata="..\\results\\rt.plt";
double gausspoint[3][2]={0,0,0,1,1,0};
double gaussconff[3]={1/3.0,1/3.0,1/3.0};


void gongetidu(double **a,double *b,double * x,int n,const double eps);


int main()
{
	
	int i,j,k,e,u,v,w,p,q;//循环变量
	int elemnum,freedom,kefreedom;
	
	double x1,y1,x2,y2,x3,y3;
	double dtemp,sum;
	ELEMENT *pe;
	Timer tm1;
	tm1.start();


	kefreedom=3;//element freedom
	elemnum=2*n*n;//
	freedom=(n+1)*(n+1);
	
	
	//网格数据
	double **pointlist;
	pointlist=new double*[(n+1)*(n+1)];
	for (i=0;i<(n+1)*(n+1);i++)
	{
		pointlist[i]=new double[2]; 
	}
	
	int **elementlist;
	elementlist=new int*[n*n];
	for (i=0;i<n*n;i++)
	{
		elementlist[i]=new int[4]; 
	}
	
	
	int **trielementlist;
	trielementlist=new int*[2*n*n];
	for (i=0;i<2*n*n;i++)
	{
		trielementlist[i]=new int[3];
		for (j=0;j<3;j++)
		{
			trielementlist[i][j]=0;
		}
	}
	
	j=0;k=0;	
	for (i=0;i<(n+1)*(n+1);i++)
	{
		pointlist[i][0]=k*(1.0/n);
		pointlist[i][1]=j*(1.0/n);
		j++;
		
		if (j>=n+1)
		{
			k++;
			j=0;
		}
	}
	j=0;k=0;
	for (i=0;i<n*n;i++)
	{
		elementlist[i][0]=j+k+n+1;
		elementlist[i][1]=j+k+n+2;	
		elementlist[i][2]=j+k+1;	
		elementlist[i][3]=j+k;
		j++;
		if (j>=n)
		{
			j=0;
			k+=n+1;
		}
	}
	for (i=0;i<2*n*n;i+=2)
	{
		trielementlist[i][0]=elementlist[i/2][0];
		trielementlist[i][1]=elementlist[i/2][3];
		trielementlist[i][2]=elementlist[i/2][2];
		
		trielementlist[i+1][0]=elementlist[i/2][0];
		trielementlist[i+1][1]=elementlist[i/2][2];
		trielementlist[i+1][2]=elementlist[i/2][1];
	}
	
	
	
	
//边界条件	
	double *zeroboundry=new double [2*n+1];
	for (i=0;i<n+1;i++)
	{
		zeroboundry[i]=i;
	}
	p=n+1;
	for (i=n+1;i<2*n+1;i++)
	{
		zeroboundry[i]=p;
		p+=(n+1);
	}
	double *nonzeroboundry=new double[2*n+1];
	double *nonzbvalue=new double [2*n+1];
	p=n;
	for (i=0;i<n+1;i++)
	{
		nonzeroboundry[i]=p;
		nonzbvalue[i]=pointlist[p][0];
		p+=n+1;
	}
	q=n*(n+1);
	for (i=n+1;i<2*n+1;i++)
	{
		nonzeroboundry[i]=q;
		nonzbvalue[i]=pointlist[q][1]*pointlist[q][1]*pointlist[q][1];
		q++;
	}
	
	
	
	//单刚	
	double **ke;
	ke=new double*[kefreedom];
	for (i=0;i<kefreedom;i++)
	{
		ke[i]=new double[kefreedom];
		for(j=0;j<kefreedom;j++)
		{
			ke[i][j]=0;
		}
	}
	
	
	//总刚
	double **KK;
	KK=new double*[freedom];
	for (i=0;i<freedom;i++)
	{
		KK[i]=new double[freedom];
		for(j=0;j<freedom;j++)
		{
			KK[i][j]=0;
		}
		
	}
	//right hand 
	double *fe=new double[kefreedom];
	double *rightf=new double[freedom];
	for (i=0;i<kefreedom;i++)
	{
		fe[i]=0;
	}
	for (i=0;i<freedom;i++)
	{
		rightf[i]=0;
	}
	tm1.stop();
	tm1.diplayelapsedtime();
	tm1.start();
	
	
//开始单元循环	
	
	for (e=0;e<elemnum;e++)
	{
		pe=new ELEMENT;
		x1=pointlist[trielementlist[e][0]][0];
		y1=pointlist[trielementlist[e][0]][1];
		x2=pointlist[trielementlist[e][1]][0];
		y2=pointlist[trielementlist[e][1]][1];
		x3=pointlist[trielementlist[e][2]][0];
		y3=pointlist[trielementlist[e][2]][1];
		
		pe->initelem(x1,y1,x2,y2,x3,y3);

		
		
		for (i=0;i<kefreedom;i++)
		{
			for (j=0;j<kefreedom;j++)
			{
				// 生成单刚矩阵
				
				sum=0;
				for (k=0;k<3;k++)
				{	
					dtemp=(pe->dshapexy(i+1,0,gausspoint[k][0],gausspoint[k][1]))*(pe->dshapexy(j+1,0,gausspoint[k][0],gausspoint[k][1]));
					dtemp+=(pe->dshapexy(i+1,1,gausspoint[k][0],gausspoint[k][1]))*(pe->dshapexy(j+1,1,gausspoint[k][0],gausspoint[k][1]));
					dtemp*=pe->getdetjcb();
					sum+=gaussconff[k]*dtemp;
				}

				ke[i][j]=0.5*sum;				
				
				//组装到总刚
				u=trielementlist[e][i];
				v=trielementlist[e][j];
				KK[u][v]+=ke[i][j];
			}
		}
		
		
		dtemp=0;
		for (i=0;i<kefreedom;i++)
		{
			sum=0;
 			for (k=0;k<3;k++)
 			{	
 				dtemp=6*(pe->zuobiaox(gausspoint[k][0],gausspoint[k][1]))*(pe->zuobiaoy(gausspoint[k][0],gausspoint[k][1]));
 				dtemp*=(pe->getdetjcb()*(pe->shape(i+1,gausspoint[k][0],gausspoint[k][1])));
 				sum+=gaussconff[k]*dtemp;
 			}

			//组装
			fe[i]=-0.5*sum;
			w=trielementlist[e][i];
			rightf[w]=fe[i];
		}
		delete pe;
	}	
	


	/*
	*	处理边界条件
	*/
	for (i=0;i<2*n+1;i++)
	{
		k=zeroboundry[i];
		for (p=0;p<freedom;p++)
		{
			KK[k][p]=0;
			KK[p][k]=0;
		}
		KK[k][k]=1;
		rightf[k]=0;
	}
	for (i=0;i<2*n+1;i++)
	{
		k=nonzeroboundry[i];
		for(p=0;p<freedom;p++)
		{
			rightf[p]-=KK[p][k]*nonzbvalue[i];
		}
		rightf[k]=nonzbvalue[i];
		for (p=0;p<freedom;p++)
		{
			KK[k][p]=0;
			KK[p][k]=0;
		}
		KK[k][k]=1;
	}


	
	
	
	//求解总刚
	double *uxy=new double[freedom];
	for (i=0;i<freedom;i++)
	{
		uxy[i]=0;
	}
	tm1.stop();
	tm1.diplayelapsedtime();
	tm1.start();
	
	gongetidu(KK,rightf,uxy,freedom,eps);				
	
	
	
//输出数据	
	ofstream rt( Resultsdata, ios_base::out );
	rt<<"TITLE=\"wangge\""<<endl
		<<"VARIABLES=x,y,u"<<endl;
	rt<<"ZONE N="<<freedom<<", E="<<elemnum<<", F=FEPOINT, ET=triangle"<<endl;
	
	for (i=0;i<freedom;i++)
	{
		rt<<pointlist[i][0]<<"\t"<<pointlist[i][1]<<"\t"<<uxy[i]<<endl;
	}
	for(i=0;i<elemnum;i++)
	{
		for (j=0;j<3;j++)
		{
			rt<<trielementlist[i][j]+1<<"\t";
		}
		rt<<endl;
	}
	double *u0xy=new double[freedom];

	rt<<"ZONE N="<<freedom<<", E="<<elemnum<<", F=FEPOINT, ET=triangle"<<endl;
	
	for (i=0;i<freedom;i++)
	{
		u0xy[i]=pow(pointlist[i][1],3)*pointlist[i][0];
		rt<<pointlist[i][0]<<"\t"<<pointlist[i][1]<<"\t"<<u0xy[i]<<endl;
	}
	for(i=0;i<elemnum;i++)
	{
		for (j=0;j<3;j++)
		{
			rt<<trielementlist[i][j]+1<<"\t";
		}
		rt<<endl;
	}
	sum=0;
	for (i=0;i<freedom;i++)
	{
		sum+=(uxy[i]-u0xy[i])*(uxy[i]-u0xy[i]);
	}
	cout<<sqrt(sum)<<endl;
	rt.close();
	tm1.stop();
	tm1.diplayelapsedtime();

	return 0;
}


void gongetidu(double **a,double *b,double * x,int n,const double eps)
{
	int i,j;
	int num;
	double alphak,betak,add,err,h,l,t;
	add=alphak=betak=0;
	
	double *p=new double [n];
	assert(p!=NULL);
	double *r=new double [n];
	assert(r!=NULL);
	
	for (i=0;i<n;i++)
	{
		p[i]=r[i]=b[i];
	}
	num=1;
	do 
	{
		//计算ALPHAK；
		h=l=0;
		for (i=0;i<n;i++)
		{
			t=0;
			for (j=0;j<n;j++)
			{
				t+=a[i][j]*p[j];
			}
			l+=p[i]*t;
			h+=r[i]*p[i];
		}
		alphak=h/l;
		/////////////////////////////////////////////		
		for (i=0;i<n;i++)
		{
			x[i]+=alphak*p[i];
		}
		
		
		for (i=0;i<n;i++)
		{
			h=0;
			for (j=0;j<n;j++)
			{
				h+=a[i][j]*x[j];
			}
			r[i]=b[i]-h;
		}
		//计算BETAK；
		h=l=0;
		for (i=0;i<n;i++)
		{	
			t=0;
			for (j=0;j<n;j++)
			{
				t+=a[i][j]*p[j];
			}
			h+=r[i]*t;
			l+=p[i]*t;
		}
		betak=-h/l;
		
		for (i=0;i<n;i++)
		{
			p[i]=r[i]+betak*p[i];
		}
		
		
		//计算误差
		err=0;
		for(i=0;i<n;i++) err+=r[i]*r[i];
		if(num%10==0)
 		{
			cout.precision(6);
			cout.setf(ios::scientific);
			cout.setf(ios::right);
			cout<<setw(5)<<setfill('#')<<num<<" "<<setw(10)<<sqrt(err)<<endl;
			cout.unsetf(ios::scientific);
		}
		num++;
	} while(sqrt(err)>eps);
	
	delete [] p;
	p=NULL;
	delete [] r;
	r=NULL;
}
