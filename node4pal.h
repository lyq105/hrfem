// node4pal.h: interface for the node4pal class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _NODE4PAL_H_
#define _NODE4PAL_H_

#include "deftype.h"
#include "AbsMtx.h"
#include "mesh.h"
#include "Point.h"


// enum Coord{CoordX,CoordY,CoordZ};
// enum RefCoord{CoordKesi,CoordYita,CoordZita};
enum WhichEage{One1,Two2,Three3,Four4};



/*!
	�Ľ�����޵�Ԫ�Ķ���	
*/

class Equation;

class node4pal  
{

public:
	node4pal();
	virtual ~node4pal();

public:
	/*!   ��ʼ����Ԫ      */
	bool Initial(element2D& ele);


	double getdetjcb(void){return detjcb;};
	double getarea(void){return area;};

	/*!	�κ����Լ��κ����ĵ��������Բο�����Ϊ�������*/

	double ShapeFunction(Point & pt,int nuNode) ;//�κ��� 
	double DShapeFunction(Point & pt,int nuNode,RefCoord nVrib);//�κ����Բο�ƽ������ĵ���
	double DShapeFunctionXY(Point & pt,int nuNode,Coord nvrib);//�κ���������ƽ������ĵ���

	/*!		����任	 */

	Point Transform(Point& pt);
	Point ReverseTransform(Point& pt);

	/*!		��Ԫ����     */

	double GaussIntOnElem(Equation& eq,int i,int j,WhichEq mark);
	FullMtx FormElemMatrix(Equation& eq,WhichEq mark);

	/*!		��Ԫ�߽����          */

	double node4pal::GaussIntOnEage(Equation& eq,int i,WhichEage mark);
	Vec node4pal::FormElemRHS(Equation& eq,WhichEq mark);

private:

	double area;
	double detjcb;
	int nNode;
	Point  Vertex[4];	
	Point  GaussPoint[4];
	double GaussWeight[4];
	FullMtx jcbdxy;
	FullMtx jcbdyk;
	FullMtx est;
	Vec		rhs;
};


#endif //_NODE4PAL_H_
