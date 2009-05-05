// FEM.h: interface for the FEM class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEM_H__DD248506_05AF_4C64_B1ED_A222A3BE342A__INCLUDED_)
#define AFX_FEM_H__DD248506_05AF_4C64_B1ED_A222A3BE342A__INCLUDED_


#include "AbsMtx.h"
#include "SparseMtxEx.h"
#include "mesh.h"
#include "node4pal.h"
#include "Equation.h"


typedef FullMtx ElementStiff;
typedef SparseMtxEx Stiffness;

typedef SparseMtxEx QualityMat;
typedef SparseMtxEx Matrix;

//typedef node4pal Element;
typedef Vec Solution;


class FEM  
{

public:
	FEM();
	virtual ~FEM();

/*
	初始化
	*/
	bool Initial(char*,char*);
/*
	形成矩阵
	*/
	bool FormMatrix(Matrix& mat,Equation & eq,WhichEq EQtype);
/*	
	组装矩阵操作
	*/
	bool AssembleMatrix(Stiffness& stf,ElementStiff& estf,Elements& ele);
/*
	右端项的相关函数
	*/
	bool FormRHS(Vec& ve,Equation & eq,WhichEq Eqtype);
    bool AssembleRHS(Vec& rhs,Vec& elevec,Elements& ele);
private:
	
	int TotleFreedom;
	Mesh Me;
};

#endif // !defined(AFX_SOLVESEQUENCE_H__DD248506_05AF_4C64_B1ED_A222A3BE342A__INCLUDED_)
