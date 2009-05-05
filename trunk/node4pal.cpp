// node4pal.cpp: implementation of the node4pal class.
//
//////////////////////////////////////////////////////////////////////

#include "node4pal.h"
#include <assert.h>
#include <math.h>



node4pal::node4pal()
//构造函数
{
	GaussWeight[0]=1;GaussWeight[1]=1;
	GaussWeight[2]=1;GaussWeight[3]=1;

	GaussPoint[0].setx(1.0/sqrt(3));GaussPoint[0].sety(1.0/sqrt(3));
	GaussPoint[1].setx(1.0/sqrt(3));GaussPoint[1].sety(-1.0/sqrt(3));
	GaussPoint[2].setx(-1.0/sqrt(3));GaussPoint[2].sety(1.0/sqrt(3));
	GaussPoint[3].setx(-1.0/sqrt(3));GaussPoint[3].sety(-1.0/sqrt(3));

}
node4pal::~node4pal()
{

}

bool node4pal::Initial(element2D& ele)
{
	double L1,L2;
	//读入单元物理坐标
	
	Vertex[0]=ele.nodes[0];
	Vertex[1]=ele.nodes[1];
	Vertex[2]=ele.nodes[2];
	Vertex[3]=ele.nodes[3];
	
	L1=0.5*(Vertex[0].getx()-Vertex[1].getx());
	L2=0.5*(Vertex[0].gety()-Vertex[3].gety());
	area=4*L1*L2;
	detjcb=fabs(L1*L2);
	
	
	//雅可比矩阵dyk/dxy
	jcbdxy.insert(0,0,1.0/L1);
	jcbdxy.insert(0,1,0);
	jcbdxy.insert(1,0,0);
	jcbdxy.insert(1,1,1.0/L2);

	//雅可比矩阵dxy/dyk
	jcbdyk.insert(0,0,L1);
	jcbdyk.insert(0,1,L1);
	jcbdyk.insert(1,0,L1);
	jcbdyk.insert(1,1,L2);

	return 1;
}


double node4pal::ShapeFunction(Point & pt,int nuNode) //形函数 
{
	switch(nuNode)
	{
	case 1:
		return 0.25*(1+pt.getx())*(1+pt.gety());
	case 2:
		return 0.25*(1-pt.getx())*(1+pt.gety());
	case 3:
		return 0.25*(1-pt.getx())*(1-pt.gety());
	case 4:
		return 0.25*(1+pt.getx())*(1-pt.gety());
	}
	return 1000;

}
double node4pal::DShapeFunction(Point & pt,int nuNode,RefCoord nVrib)//形函数对参考平面坐标的导数
{
	switch(nuNode)
	{	
	case 1:
		if (nVrib==CoordKesi) 
			{return 0.25*(1+pt.gety());}
		else return 0.25*(1+pt.getx());
	case 2:
		if(nVrib==CoordKesi)
			{return -0.25*(1+pt.gety());}
		else return 0.25*(1-pt.getx());
	case 3:
		if(nVrib==CoordKesi)
			{return -0.25*(1-pt.gety());}
		else return -0.25*(1-pt.getx());
	case 4:
		if(nVrib==CoordKesi)
			{return 0.25*(1-pt.gety());}
		else return -0.25*(1+pt.getx());
	}
	return 1000;
}
double node4pal::DShapeFunctionXY(Point & pt,int nuNode,Coord nvrib)//形函数对物理平面坐标的导数
{
	double L1,L2;

	L1=0.5*(Vertex[0].getx()-Vertex[1].getx());
	L2=0.5*(Vertex[0].gety()-Vertex[3].gety());
	
	if(nvrib==CoordX)return 1.0/L1*DShapeFunction(pt,nuNode,CoordKesi);
	else return 1.0/L2*DShapeFunction(pt,nuNode,CoordYita);
}

/*!
		物理坐标到参考坐标	
  */

Point node4pal::Transform(Point& pt)
{
	return Point(0,0,0);
}

/*!
		参考坐标到物理坐标	
		(x,y):=RT(kesi,yita)
  */

Point node4pal::ReverseTransform(Point& refpt)
{
	Point p;
	double L1,L2,x0,y0;
	L1=0.5*(Vertex[0].getx()-Vertex[1].getx());
	L2=0.5*(Vertex[0].gety()-Vertex[3].gety());
	x0=0.5*(Vertex[0].getx()+Vertex[1].getx());
	y0=0.5*(Vertex[0].gety()+Vertex[3].gety());
	
	p.setx(L1*refpt.getx()+x0);
	p.sety(L2*refpt.gety()+y0);
	
	return  p;

}

double node4pal::GaussIntOnElem(Equation& eq,int i,int j,WhichEq mark)
{
	//高斯积分
	double sum=0;
	for (int p=0;p<4;p++)sum+=GaussWeight[i]*eq(*this,GaussPoint[i],i,j);
	return getdetjcb()*sum; 
}


FullMtx node4pal::FormElemMatrix(Equation& eq,WhichEq mark)
{
	FullMtx stiff(4,4);
	for (int i=0;i<4;i++)
	{
		for (int j=0;j<4;j++)
		{
			stiff.insert(i,j,GaussIntOnElem(eq,i,j,mark));
		}
	}
	return stiff; 
}

double node4pal::GaussIntOnEage(Equation& eq,int i,WhichEage mark)
{
	//高斯积分
	double sum=0;
	for (int p=0;p<4;p++)sum+=GaussWeight[i]/*eq(*this,GaussPoint[i],i)*/;
	return getdetjcb()*sum; 
}

Vec node4pal::FormElemRHS(Equation& eq,WhichEq mark)
{
	double EleVec[4];
	
	for (int i=0;i<4;i++)
	{
		EleVec[i]=GaussIntOnEage(eq,i,mark);	
	}

	return Vec(4,EleVec);
}

