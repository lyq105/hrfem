// FEM.cpp: implementation of the FEM class.
//
//////////////////////////////////////////////////////////////////////

#include "SolveSequence.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEM::FEM()
{
}

FEM::~FEM()
{
	Me.ReleaseMesh();
}

//初始化网格等计算数据
bool FEM::Initial(char* ElementFilename,char* NodesFileName)
{

	Me.InitialMesh(ElementFilename,NodesFileName);
	TotleFreedom=Me.nNode;
	return true;	
};


bool FEM::FormMatrix(Matrix& mat,Equation & eq,WhichEq EQtype)
{
	for (int i=0;i<Me.nElem;i++)
	{
		Element* pEt=new Element;
		pEt->Initial(Me.pElements[i]);
		AssembleMatrix(mat,pEt->FormElemMatrix(eq,EQtype),(*Me.pElements));
		delete pEt;
	}
	return true;
};


bool FEM::AssembleMatrix(Stiffness& stf,ElementStiff& estf,Elements& ele)
{
	for (int i=0;i<ele.nNodes;i++)
	{
		for (int j=0;j<ele.nNodes;j++)
		{
			double val=stf(ele.i,ele.j)+estf[i][j];
			stf.insert(ele.i,ele.j,val);
		}
	}
	return true;
}
/*----------------------------------------------------*/

/*-------------------------右端项---------------------*/

bool FEM::FormRHS(Vec& ve,Equation & eq,WhichEq Eqtype)
{
	for (int i=0;i<Me.nElem;i++)
	{
		Element* pEt=new Element;
		pEt->Initial(Me.pElements[i]);
		AssembleRHS(ve,pEt->FormElemRHS(eq,Eqtype),(*Me.pElements));
		delete pEt;
	}
	return true;
};

/*************************/

bool FEM::AssembleRHS(Vec& rhs,Vec& elevec,Elements& ele)
{
	for (int i=0;i<ele.nNodes;i++)
	{
			rhs[ele.i]+=elevec[i];
	}
	return true;
}

