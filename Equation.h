// Equation.h: interface for the Equation class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_EQUATION_H__53CB1C24_E448_48A4_9294_8840D511A11A__INCLUDED_)
#define AFX_EQUATION_H__53CB1C24_E448_48A4_9294_8840D511A11A__INCLUDED_


#include "AbsMtx.h"
#include "Point.h"
#include "deftype.h"
#include "defElement.h"

class Equation  
{
public:
	Equation();
	virtual ~Equation();
	
	double operator()(Point&,Element&,WhichEq);

	double functional1(Point&,Element&,int i,int j);
	double functional2(Point&,Element&,int i);

};

#endif // !defined(AFX_EQUATION_H__53CB1C24_E448_48A4_9294_8840D511A11A__INCLUDED_)
