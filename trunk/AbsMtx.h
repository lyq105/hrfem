/**************************************************************************
*
*				Linear Algebra's Definition
*
*					2008.9.12
*					    0.2
*		        Modified 2009.3.24
****************************************************************************/




#ifndef _MTX_H_
#define _MTX_H_	
#include <iostream>
using std::ostream;
using std::istream;


/*--------------------------------------------------------
|
|				Vector Class'Definition
|
---------------------------------------------------------*/ 
class Vec
{
public:
//construct function
	
	Vec(int,const double *);
	Vec(int = 1,double =0);
	Vec(const Vec&);

//Destruct function

	~Vec(){delete [] vr;}

//Interface of Vec Class
public:

	Vec &operator=(const Vec&);	
	Vec &operator+=(const Vec&);
	Vec &operator-=(const Vec&);
	double &operator [] (int i){return vr[i];}
	double &operator [] (int i)const{return vr[i];}

//
public:

	double maxnorm() const;
	double twonorm() const;
	int size()const {return lenth;}	
	bool Setlenth(int len);
//friend function
	friend double dot(const Vec &,const Vec &);
	friend Vec operator+(const Vec&);
	friend Vec operator-(const Vec&);
	friend Vec operator+(const Vec&,const Vec&);
	friend Vec operator-(const Vec&,const Vec&);
	friend Vec operator*(double,const Vec&);
	friend Vec operator*(const Vec&,double);
	friend Vec operator*(const Vec&,const Vec&);
	friend Vec operator/(const Vec&,double);
	friend ostream& operator << (ostream&,const Vec&);

private:
	int lenth;		 //vector's lenth
	double * vr;     //values of vector
};
/*--------------------------------------------------------------
|
|				Abstract Matrix' Definition
|
--------------------------------------------------------------*/ 
class AbsMtx  
{
public:
	AbsMtx(){};
	virtual ~AbsMtx(){};
public:
	virtual Vec operator*(const Vec &) const = 0;
//conjugate 
	int CG(Vec &,const Vec&,double &,int &,int);
// 	int GAUSSELIMATE();
//      int GMRERS();
protected:
	int nrows;
	virtual Vec preconding(const Vec & r,int i=0) const = 0;
};

/*-----------------------------------------------------------------
|
|				  Dense Matrix 'Definition
|
-----------------------------------------------------------------*/
class FullMtx:public AbsMtx
{
public:
//construct	
	FullMtx(int n,int m,double**);
	FullMtx(int n=1,int m=1,double t = 0);

	FullMtx(const FullMtx &);

//destruct
	~FullMtx()
	{
		for (int i=0;i<nrows;i++)
		{
			delete[] mx[i];
		}
		delete mx;
	}
//interface 
public:
	friend ostream & operator << (ostream &,FullMtx &);
// overloaded operators
	FullMtx & operator=(const FullMtx &);
	FullMtx & operator+=(const FullMtx &);
	FullMtx & operator-=(const FullMtx &);
	FullMtx  operator +(const FullMtx &);
	FullMtx  operator -(const FullMtx &);
	FullMtx & operator+();
	FullMtx  operator-();
	inline double* operator[](int i){return mx[i];}
	inline double* operator[](int i)const{return mx[i];}
	Vec operator*(const Vec&) const;
	Vec operator*(const Vec&);
// Gauss Elimination Method
	void GaussElim(Vec &) const;
// row and column
	inline int row(){return nrows;}
	inline int column(){return ncols;}
	int insert(int i,int j,double val);
	int Setsize(int n,int m);
private:
	int ncols;
	double ** mx;
	Vec preconding(const Vec & r,int i=0) const;
};

/*-----------------------------------------------------------------
|
|				Sparse Matrix Class'Definition
|
-----------------------------------------------------------------*/
class SparseMtx:public AbsMtx
{
public:
	friend ostream & operator << (ostream &,SparseMtx &);

	//construction and destruction 
	SparseMtx(int nrow,int len,double *t,int *c,int *f);
	SparseMtx(int nrow,int len);
	SparseMtx(const SparseMtx &);
	SparseMtx(FullMtx & Fmt);

	~SparseMtx(){delete [] sra;delete [] fnz;delete [] clm;}

public:
	//overload operator
	SparseMtx & operator=(const SparseMtx &);
	Vec operator*(const Vec&) const;
	double operator[](int i)const{return sra[i];}
	double operator()(int i,int j);
	//other
	int getclm(int i) const {return clm[i];}
	int getfnz(int i) const {return fnz[i];}
	double getsra(int i) const {return sra[i];}

	int row()const{return nrows;}
	int getlenth() const {return lenth;}
	int getlenth() {return lenth;}

private:
	int lenth;
	int ncols;
	int * clm;
	int * fnz;
	double * sra;
	Vec preconding(const Vec &,int i = 0) const;

};

#endif //_MTX_H_
