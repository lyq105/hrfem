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
	四结点有限单元的定义	
*/

class Equation;

class node4pal  
{

public:
	node4pal();
	virtual ~node4pal();

public:
	/*!   初始化单元      */
	bool Initial(element2D& ele);


	double getdetjcb(void){return detjcb;};
	double getarea(void){return area;};

	/*!	形函数以及形函数的导数，均以参考坐标为输入参数*/

	double ShapeFunction(Point & pt,int nuNode) ;//形函数 
	double DShapeFunction(Point & pt,int nuNode,RefCoord nVrib);//形函数对参考平面坐标的导数
	double DShapeFunctionXY(Point & pt,int nuNode,Coord nvrib);//形函数对物理平面坐标的导数

	/*!		坐标变换	 */

	Point Transform(Point& pt);
	Point ReverseTransform(Point& pt);

	/*!		单元积分     */

	double GaussIntOnElem(Equation& eq,int i,int j,WhichEq mark);
	FullMtx FormElemMatrix(Equation& eq,WhichEq mark);

	/*!		单元边界积分          */

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
