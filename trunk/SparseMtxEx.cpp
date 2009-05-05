#include "SparseMtxEx.h"
#include <assert.h>
#include <math.h>
/*-------------------------------------------------------
|
|			Sparse Matrix Class
|
--------------------------------------------------------*/

SparseMtxEx::SparseMtxEx(int rowsnum,int len,double *t,int *c,int *f)
{
	int i;
	nrows = rowsnum;
	lenth = len;
	fnz=new int[nrows+1];

	for ( i=0;i<lenth;i++)
	{
		sra.push_back(t[i]);
		clm.push_back(c[i]);
	}
	for( i=0;i<nrows;i++)fnz[i]=f[i];
	fnz[nrows]=lenth;
}

SparseMtxEx::SparseMtxEx(int rows)
{
	int i;
	nrows = rows;
	lenth = rows;
	fnz=new int[nrows+1];

	for ( i=0;i<lenth;i++)
	{
		sra.push_back(1);
		clm.push_back(i);
	}
	for( i=0;i<nrows;i++)fnz[i]=i;
	fnz[nrows]=lenth;


	// 	for (int k=0;k<lenth;k++)
	// 	{
	// 		cout<<k<<"---"<<sra[k]<<"___"<<clm[k]<<"\n";
	// 	}
	// 
	// 	for (int l=0;l<nrows;l++)
	// 	{
	// 		cout<<fnz[l]<<"\n";
	// 	}

}
SparseMtxEx::SparseMtxEx(const SparseMtxEx &spm)
{
	int i;
	nrows = spm.nrows;
	lenth = spm.lenth;
	fnz=new int[nrows+1];

	for ( i=0;i<lenth;i++)
	{
		sra.push_back(spm.sra[i]);
		clm.push_back(spm.clm[i]);
	}
	for( i=0;i<=nrows;i++)fnz[i]=spm.fnz[i];
}

SparseMtxEx::SparseMtxEx(SparseMtx &spm)
{
	int i;
	nrows = spm.row();
	lenth = spm.getlenth();
	fnz=new int[nrows+1];

	for (i=0;i<lenth;i++)
	{
		sra.push_back(spm.getsra(i));
		clm.push_back(spm.getclm(i));
	}
	for(i=0;i<=nrows;i++)fnz[i]=spm.getfnz(i);
}

SparseMtxEx::SparseMtxEx(FullMtx &fmt)
{
	int  i,j;
	int nNonzero=0;
	for (i=0;i<fmt.row();i++)
		for(j=0;j<fmt.row();j++)
			if(fmt[i][j]!=0)nNonzero++;

	nrows=fmt.row();
	lenth=nNonzero;

	fnz=new int[nrows+1];


	int p=0;
	for ( i=0;i<fmt.row();i++)
	{	
		int k=0;
		for(j=0;j<fmt.row();j++)
			if(fmt[i][j]!=0)
			{
				if(k==0){fnz[i]=p;k=1;}
				sra.push_back(fmt[i][j]);
				clm.push_back(j);
			}
	}
	fnz[nrows]=lenth;
}

Vec SparseMtxEx::operator * (const Vec &v)const
{
	assert(nrows==v.size());
	Vec tm(nrows);
	for (int i=0;i<nrows;i++)
	{
		for (int j=fnz[i];j<fnz[i+1];j++)
		{
			tm[i]+=sra[j]*v[clm[j]];
		}
	}
	return tm;
}

SparseMtxEx& SparseMtxEx::operator = (const SparseMtxEx& spm)
{
	int i;
	if(this!=&spm)
	{
		assert(nrows==spm.nrows);	
		for (i=0;i<spm.getlenth();i++)
		{
			sra.push_back(spm.sra[i]);
			clm.push_back(spm.clm[i]);
		}
		for(i=0;i<=nrows;i++)fnz[i]=spm.fnz[i];
	}
	return *this;
}

double SparseMtxEx::operator ()(int i,int j)
{
	for(int p=fnz[i];p<fnz[i+1];p++)
		if (clm[p]==j) return sra[p];
	return 0;
}

Vec SparseMtxEx::preconding(const Vec &r,int precn /* = 0 */)const
{
	int i,j;
	if (precn==0)
	{
		return r;
	}
	else if (precn ==1)
	{
		Vec diag(nrows);
		for(i=0;i<nrows;i++)
		{
			for(j =fnz[i];j<fnz[i+1];j++)
			{
				diag[i]+=sra[j]*sra[j];
			}
			diag[i]=sqrt(diag[i]);
		}
		Vec z(nrows);
		for (i=0;i<nrows;i++)z[i]=r[i]/diag[i];
		return z;
	}else if(precn==2)
	{
		const double omega =1.2;

		Vec diag(nrows);
		for (i=0;i<nrows;i++)
		{
			for (j = fnz[i];j<fnz[i+1];j++)
			{
				if(clm[j]==i)diag[i]+=sra[j]*sra[j];
				else diag[i]+=omega*omega*sra[i]*sra[j];
			}
			diag[i] = sqrt(diag[i]);
		}

		Vec z(nrows);
		for (i=0;i<nrows;i++)
		{
			double sum = 0;
			for(j =fnz[i];j<fnz[i+1];j++)
				if (clm[j]<i)
				{
					sum+=sra[j]*z[clm[j]];
				}
				z[i]=(r[i]-omega*sum)/diag[i];
		}
		for (i=nrows-1;i>=0;i--)
		{
			double sum=0;
			for(j=fnz[i];j<fnz[i+1];j++)
			{
				if(clm[j]>i)sum+=sra[j]*z[clm[j]];
			}
			z[i]-=omega*sum/diag[i];
		}
		return z;
	}else{exit(1);}
}
ostream & operator << (ostream & ostm,SparseMtxEx & spt)
{
	int i,j;
	for (i =0; i<spt.row();i++)
	{
		for (j =0; j<spt.row();j++)
		{
			ostm<<spt(i,j)<<"\t";
		}
		ostm<<"\n";
	}
	return ostm;
}

int SparseMtxEx::insert(int i,int j,double value)
{
	assert(i<nrows&&j<nrows);
	if (value==0)//插入的元素是0
		{
			if((*this)(i,j)==0)return 1;
			else//插入位置上是非零元
			{
				for(int p=fnz[i];p<fnz[i+1];p++)
					if (clm[p]==j) 
					{
						sra.erase(sra.begin()+p);
						clm.erase(clm.begin()+p);
						break;
					}
					//更新标号
				for (int q=i+1;q<=nrows;q++)
					{
						fnz[q]-=1;
					}

			}
		}

	else //插入的元素不是0
	{
		if ((*this)(i,j)!=0)// 插入位置上的元素是非零的
		{
			for(int p=fnz[i];p<fnz[i+1];p++)
				if (clm[p]==j) {sra[p]=value;break;}
			return 1;
		}

		else				// 插入位置上的元素是零
		{
			if (j<clm[fnz[i]])// 在第一个非零元之前？
			{
				sra.insert(sra.begin()+fnz[i],value);
				clm.insert(clm.begin()+fnz[i],j);
			}
			else{
				if (j>clm[fnz[i+1]-1])//在最后一个非零元之后？
				{
					sra.insert(sra.begin()+fnz[i+1]-1,value);
					clm.insert(clm.begin()+fnz[i+1]-1,j);
				}

				else//中间位置
				{
				for(int p=fnz[i];p<fnz[i+1];p++)
				{
					if (j<clm[p])
					{
						sra.insert(sra.begin()+p-1,value);
						clm.insert(clm.begin()+p-1,j);
						break;
					}
				}
				}
			}
			//更新非零元标号
			for (int q=i+1;q<=nrows;q++)
			{
				fnz[q]+=1;
			}
		}
	}
	return 1;
}