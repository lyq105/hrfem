/**************************************************************
*
*             
*
*This is The linear Algra 's implatement			
*
*
*
*2008.11.26 modified
*
***************************************************************/
#include "AbsMtx.h"
#include <assert.h>
#include <cmath>
#include <iostream>
using std::cout;
using std::endl;


/*-------------------------------------------------------------
|        
|				 Vector class
|
--------------------------------------------------------------*/


//construction function
Vec::Vec(int n,const double * abd)
{
	lenth = n;
	vr = new double[lenth];
	for (int i=0;i<lenth;i++){vr[i]=abd[i];}
}
Vec::Vec(const Vec& vf)
{
	lenth=vf.lenth;
	vr = new double[lenth];
	for (int i=0;i<lenth;i++){vr[i]=vf.vr[i];}
}

Vec::Vec(int n,const double valu)
{
	lenth=n;
	vr = new double[lenth];
	for (int i=0;i<lenth;i++){vr[i]=valu;}
}

//operator overload
Vec& Vec::operator =(const Vec& v)
{
	if (this != &v)
	{
		assert(lenth==v.lenth);
		for(int i=0;i<lenth;i++){vr[i]=v.vr[i];}
	}
	return *this;
}
Vec & Vec::operator+=(const Vec& v)
{
	assert(lenth==v.lenth);
	for(int i=0;i<lenth;i++){vr[i]+=v.vr[i];}
	return *this;
}
Vec & Vec::operator-=(const Vec& v)
{
	assert(lenth==v.lenth);
	for(int i=0;i<lenth;i++){vr[i]-=v.vr[i];}
	return *this;
}

//norm
double Vec::maxnorm()const
{
	double nm = fabs(vr[0]);
	for (int i = 1;i<lenth;i++)
	{
		if (nm<=fabs(vr[i]))
		{
			nm=fabs(vr[i]);
		}
	}
	return nm;
}
double Vec::twonorm()const
{
	double norm2 = vr[0]*vr[0];
	for (int i = 1;i<lenth;i++){norm2+=vr[i]*vr[i];}
	return sqrt(norm2);
}

bool Vec::Setlenth(int len)
{
	lenth=len;
	vr = new double[lenth];
	for (int i=0;i<lenth;i++){vr[i]=0;}
	return true;
}



// Vector's Friend Functions

Vec operator+(const Vec & v){return v;}
Vec operator-(const Vec & v){return Vec(v.lenth)-v;}

Vec operator+(const Vec & v1,const Vec & v2)
{
	assert(v1.lenth==v2.lenth);
	Vec sum = v1;
	sum += v2;
	return sum;
}
Vec operator-(const Vec & v1,const Vec & v2)
{
	assert(v1.lenth==v2.lenth);
	Vec ms = v1;
	ms -= v2;
	return ms;
}
Vec operator*(double scalar, const Vec & v)
{
	Vec tm(v.lenth);
	for (int i=0;i<v.lenth;i++)
	{
		tm[i]=scalar*v[i];
	}
	return tm;
}
Vec operator*(const Vec & v,double scalar)
{
	return scalar*v;
}

Vec operator*(const Vec&v1,const Vec&v2)
{
	assert(v1.lenth==v2.lenth);
	Vec tm(v1.lenth);
	for (int i=0;i<v1.lenth;i++)
	{
		tm[i]=v1[i]*v2[i];
	}
	return tm;
}
Vec operator / (const Vec &v,double scalar)
{
	assert(scalar!=0);
	return (1.0/scalar)*v;
}
ostream& operator << (ostream& ostr,const Vec& vt)
{
	for (int i=0;i<vt.lenth;i++)
	{
		ostr<<vt[i]<<"\t";
		if (i%10 == 5){ostr<<"\n";}
	}
	return ostr;
}

double dot(const Vec & v1,const Vec & v2)
{
	assert(v1.lenth == v2.lenth);
	double tm = v1[0]*v2[0];
	for (int i = 1;i<v1.lenth;i++)
	{
		tm+=v1[i]*v2[i];
	}
	return tm;
}

/*-------------------------------------------------------------
|        
|				 Abstract Matrix 
|
--------------------------------------------------------------*/

int AbsMtx::CG(Vec & x,const Vec & b,double & eps,int & iter,int pcn)
{
	assert(nrows == b.size());
	const int maxiter = iter;
	Vec r = b - (*this)*x;
	Vec z = preconding(r,pcn);
	Vec p = z;
	double zr = dot(z,r);
	const double stp = eps * b.twonorm();
	
	if (!r.maxnorm())
	{
		eps = 0.0;
		iter = 0;
		return 0;
	}
	
	for(iter = 0;iter<maxiter;iter++)
	{
		Vec mp =(*this)*p;
		double pap=dot(mp,p);
		double alpha = zr/pap;
		x += alpha*p;
		r -= alpha*mp;
		if (r.twonorm() <= stp)break;
		z = preconding(r,pcn);
		double zrold = zr;
		zr = dot(z,r);
		double beta = zr/zrold;
		p = z+beta*p;
	}
	
	eps = r.twonorm();
	if (iter == maxiter)return 1;
	else return 0;
}

/*-------------------------------------------------------------
|        
|				Full matrix 
|
--------------------------------------------------------------*/
FullMtx::FullMtx(int n,int m,double** dbp)
{
	nrows = n;
	ncols = m;
	mx = new double*[nrows];
	for (int i=0;i<nrows;i++)
	{
		mx[i]=new double[ncols];
		for(int j=0;j<ncols;j++)mx[i][j]=dbp[i][j];
	}
}

FullMtx::FullMtx(int n,int m,double t /* = 0 */)
{
	nrows = n;
	ncols = m;
	mx = new double*[nrows];
	for (int i=0;i<nrows;i++)
	{
		mx[i]=new double[ncols];
		for(int j=0;j<ncols;j++)mx[i][j]=t;
	}
}
FullMtx::FullMtx(const FullMtx &fmt)
{
	nrows = fmt.nrows;
	ncols = fmt.ncols;
	mx = new double*[nrows];
	for (int i=0;i<nrows;i++)
	{
		mx[i] = new double [ncols];
		for (int j = 0;j<ncols;j++){mx[i][j]=fmt.mx[i][j];}
	}
}

FullMtx& FullMtx::operator = (const FullMtx& fmt)
{
	if (this != & fmt)
	{
		assert(nrows==fmt.nrows&&ncols==fmt.ncols);
		for (int i = 0;i < nrows;i++)
			for(int j = 0;j < ncols; j++)
				mx[i][j]=fmt.mx[i][j];
	}
	return *this;
}

FullMtx& FullMtx::operator += (const FullMtx& fmt)
{
	assert(nrows==fmt.nrows&&ncols==fmt.ncols);
	for (int i = 0;i < nrows;i++)
		for(int j = 0;j < ncols; j++)
			mx[i][j]+=fmt.mx[i][j];
		return *this;
}

FullMtx& FullMtx::operator -=(const FullMtx& fmt)
{
	assert(nrows==fmt.nrows&&ncols==fmt.ncols);
	for (int i = 0;i < nrows;i++)
		for(int j = 0;j < ncols; j++)
			mx[i][j]-=fmt.mx[i][j];
		return *this;
}
FullMtx& FullMtx::operator +()
{
	return *this;
}

FullMtx FullMtx::operator -()
{
	FullMtx zero(nrows,ncols);
	return zero - *this;
}

FullMtx FullMtx::operator +(const FullMtx & fmt)
{
 	FullMtx sum =*this;
	sum +=fmt;
 	return sum;
}
FullMtx FullMtx::operator -(const FullMtx & fmt)
{
	FullMtx sum =*this;
 	sum -=fmt;
 	return sum;
}
Vec FullMtx::operator*(const Vec& vt) 
{
	assert(ncols==vt.size());
	Vec tm(nrows);
	for (int i=0;i<nrows;i++)
		for(int j=0;j<ncols;j++)
			tm[i]+=mx[i][j]*vt[j];
	return tm;
}
Vec FullMtx::operator*(const Vec& vt) const
{
	assert(ncols==vt.size());
	Vec tm(nrows);
	for (int i=0;i<nrows;i++)
		for(int j=0;j<ncols;j++)
			tm[i]+=mx[i][j]*vt[j];
		return tm;
}


Vec FullMtx::preconding(const Vec& r,int precn)const
{
	int i=0;
	if (precn==0)
	{
		return r;
	}
	else if(precn ==1)
	{
		Vec z(nrows);
		for ( i = 0; i<nrows;i++)
		{
			z[i]=r[i]/mx[i][i];
		}
		return z;
	}
	else if (precn ==2)
	{
		const double omega=1.2;
		Vec z(nrows);
		for ( i = 0;i<nrows;i++)
		{
			double sum = 0;
			for(int j = 0;j<i;j++)sum+=mx[i][j]*z[j];
			z[i] = (r[i]-omega*sum)/mx[i][i];
		}
		for( i = nrows-1;i>=0;i--)
		{
			double sum=0;
			for (int j = i+1;j<nrows;j++)sum+=mx[i][j]*z[j];
			z[i]-= omega*sum/mx[i][i];
		}
		return z;
	}else {exit(1);}
}
// Friend Methods
ostream & operator << (ostream & ostm,FullMtx & fmt)
{
	for (int i =0; i<fmt.row();i++)
	{
		for (int j =0; j<fmt.column();j++)
		{
			ostm<<fmt[i][j]<<"\t";
		}
		ostm<<"\n";
	}
	return ostm;
}
int FullMtx::insert(int i,int j,double val)
{
	assert(i<nrows&&j<ncols);
	mx[i][j]=val;
	return 1;
}
int FullMtx::Setsize(int n,int m)
{
	nrows = n;
	ncols = m;
	mx = new double*[nrows];
	for (int i=0;i<nrows;i++)
	{
		mx[i]=new double[ncols];
		for(int j=0;j<ncols;j++)mx[i][j]=0;
	}
return 0;
}

/*-------------------------------------------------------
|
|			Sparse Matrix Class
|
--------------------------------------------------------*/

SparseMtx::SparseMtx(int rowsnum,int len,double *t,int *c,int *f)
{
	int i=0;
	nrows = rowsnum;
	lenth = len;
	sra = new double [lenth];
	clm = new int [lenth];
	fnz = new int [nrows+1];
	
	for ( i=0;i<lenth;i++)
	{
		sra[i]=t[i];
		clm[i]=c[i];
	}
	for( i=0;i<nrows;i++)fnz[i]=f[i];
	fnz[nrows]=lenth;
}

SparseMtx::SparseMtx(int rows,int len)
{
	int i=0;
	nrows = rows;
	lenth = len;
	sra = new double [lenth];
	clm = new int [lenth];
	fnz = new int [nrows+1];
	
	for ( i=0;i<lenth;i++)
	{
		sra[i]=0;
		clm[i]=0;
	}
	for(  i=0;i<=nrows;i++)fnz[i]=0;
}
SparseMtx::SparseMtx(const SparseMtx &spm)
{
	int i=0;
	nrows = spm.nrows;
	lenth = spm.lenth;
	sra = new double [lenth];
	clm = new int [lenth];
	fnz = new int [nrows+1];
	
	for ( i=0;i<lenth;i++)
	{
		sra[i]=spm.sra[i];
		clm[i]=spm.clm[i];
	}
	for( i=0;i<=nrows;i++)fnz[i]=spm.fnz[i];
}

SparseMtx::SparseMtx(FullMtx &fmt)
{
	int i,j;
	int nNonzero=0;
	for ( i=0;i<fmt.row();i++)
		for( j=0;j<fmt.row();j++)
			if(fmt[i][j]!=0)nNonzero++;
	nrows=fmt.row();
	lenth=nNonzero;
	
	sra = new double [lenth];
	clm = new int [lenth];
	fnz = new int [nrows+1];
	int p=0;
	for (i=0;i<fmt.row();i++)
	{	
		int k=0;
		for(j=0;j<fmt.row();j++)
			if(fmt[i][j]!=0)
			{
				if(k==0){fnz[i]=p;k=1;}
				sra[p]=fmt[i][j];
				clm[p]=j;
				p++;
			}
	}
	fnz[nrows]=lenth;
}
Vec SparseMtx::operator * (const Vec &v)const
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
SparseMtx& SparseMtx::operator = (const SparseMtx& spm)
{
	int i=0;
	if(this!=&spm)
	{
		assert(nrows==spm.nrows&&lenth==spm.lenth);	
		for ( i=0;i<lenth;i++)
		{
			sra[i]=spm.sra[i];
			clm[i]=spm.clm[i];
		}
		for( i=0;i<=nrows;i++)fnz[i]=spm.fnz[i];
	}
	return *this;
}

double SparseMtx::operator ()(int i,int j)
{
	for(int p=fnz[i];p<fnz[i+1];p++)
		if (clm[p]==j) return sra[p];
	return 0;
}
Vec SparseMtx::preconding(const Vec &r,int precn /* = 0 */)const
{
	int i,j;
	if (precn==0)
	{
		return r;
	}
	else if (precn ==1)
	{
		Vec diag(nrows);
		for( i=0;i<nrows;i++)
		{
			for( j =fnz[i];j<fnz[i+1];j++)
			{
				diag[i]+=sra[j]*sra[j];
			}
			diag[i]=sqrt(diag[i]);
		}
		Vec z(nrows);
		for ( i=0;i<nrows;i++)z[i]=r[i]/diag[i];
		return z;
	}else if(precn==2)
	{
		const double omega =1.2;

		Vec diag(nrows);
		for ( i=0;i<nrows;i++)
		{
			for ( j = fnz[i];j<fnz[i+1];j++)
			{
				if(clm[j]==i)diag[i]+=sra[j]*sra[j];
				else diag[i]+=omega*omega*sra[i]*sra[j];
			}
			diag[i] = sqrt(diag[i]);
		}

		Vec z(nrows);
		for ( i=0;i<nrows;i++)
		{
			double sum = 0;
			for( j =fnz[i];j<fnz[i+1];j++)
				if (clm[j]<i)
				{
					sum+=sra[j]*z[clm[j]];
				}
				z[i]=(r[i]-omega*sum)/diag[i];
		}
		for ( i=nrows-1;i>=0;i--)
		{
			double sum=0;
			for( j=fnz[i];j<fnz[i+1];j++)
			{
				if(clm[j]>i)sum+=sra[j]*z[clm[j]];
			}
			z[i]-=omega*sum/diag[i];
		}
		return z;
	}else{exit(1);}
}
ostream & operator << (ostream & ostm,SparseMtx & spt)
{
	for (int i =0; i<spt.row();i++)
	{
		for (int j =0; j<spt.row();j++)
		{
			ostm<<spt(i,j)<<"\t";
		}
		ostm<<"\n";
	}
	return ostm;
}
