#ifndef _SPARSEMTXEX_H_
#define _SPARSEMTXEX_H_

#include "AbsMtx.h"

#include <vector>
using std::vector;


//ϡ�跽��
class SparseMtxEx:public AbsMtx 
{
public:
	friend ostream & operator << (ostream &,SparseMtxEx &);
	//construction and destruction 
	SparseMtxEx(int nrow,int len,double *t,int *c,int *f);
	SparseMtxEx(int nrow = 1);
	SparseMtxEx(const SparseMtxEx &);
	SparseMtxEx(SparseMtx &);
	SparseMtxEx(FullMtx & Fmt);
	//
	~SparseMtxEx(){delete [] fnz;}
public:
	//overload operator
	SparseMtxEx & operator=(const SparseMtxEx &);
	Vec operator*(const Vec&) const;
	double operator[](int i)const{return sra[i];}
	double operator()(int i,int j);

	bool Setrow(int nrow){nrows=nrow;return true;}

	//other
	int getclm(int i) const {return clm[i];}
	int getfnz(int i) const {return fnz[i];}
	double getsra(int i) const {return sra[i];}

	int row()const{return nrows;}
	int getlenth(){return lenth;}
	int getlenth() const {return lenth;}

	int insert(int i,int j,double value);

private:
	int lenth;
	int* fnz;			//ÿ�е�һ������Ԫ����sra�е�λ��
	vector<int> clm;    //sra����Ӧλ�õ�Ԫ�ص��к�
	vector<double> sra; //����Ԫ�ص�ֵ
	Vec preconding(const Vec &,int i = 0) const;
};


#endif //_SPARSEMTXEX_H_