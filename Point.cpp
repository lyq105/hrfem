/////////////////////////////////////////////////////////////////
//Class Point
/////////////////////////////////////////////////////////////////



#include "Point.h"


Point::Point(double a/* =0 */,double b/* =0 */,double c/* =0 */)
{
	x=a;
	y=b;
	z=c;
}
Point::Point(Point & pt)
{
	x=pt.getx();
	y=pt.gety();
	z=pt.getz();
}
Point& Point::operator = (const Point &pt)
{
	if(this != & pt)
	{
		x=pt.getx();
		y=pt.gety();
		z=pt.getz();
	}
	return *this;
}

Point operator *(Point & pt ,double rhs)
{
	Point ps;
	ps.setx(pt.getx()*rhs);
	ps.sety(pt.gety()*rhs);
	ps.setz(pt.getz()*rhs);
	return ps;
}
Point operator *(double lhs,Point & pt )
{
	return pt*lhs;
}

Point& Point::operator +=(const Point &pt)
{
	x+=pt.getx();
	y+=pt.gety();
	z+=pt.getz();
	return *this;
}

Point operator+(const Point & lhs,const Point &rhs)
{
	Point spt;
	spt.setx(lhs.x+rhs.x);
	spt.sety(lhs.y+rhs.y);
	spt.setz(lhs.z+rhs.z);
	return spt;
}

Point operator-(const Point & lhs,const Point& rhs)
{
	Point spt;
	spt.setx(lhs.x-rhs.x);
	spt.sety(lhs.y-rhs.y);
	spt.setz(lhs.z-rhs.z);
	return spt;
}