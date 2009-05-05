#ifndef _POINT_H_
#define _POINT_H_

class Point
{
public:
	Point(Point&);
	Point(double =0,double =0,double =0);
	~Point(){};
public:
	//	返回坐标值
	inline double getx() {return x;}
	inline double gety() {return y;}
	inline double getz() {return z;}
	//  重载const函数
	inline double getx()const{return x;}
	inline double gety()const{return y;}
	inline double getz()const{return z;}
	
	//	确定坐标值
	inline void setx(double a){ x = a;}
	inline void sety(double b){ y = b;}
	inline void setz(double c){ z = c;}
	//	运算符重载
	Point &operator = (const Point & pt);
	Point &operator +=(const Point & pt);
	
	friend Point operator * (Point& pt,double rhs);
	friend Point operator * (double lhs,Point &pt);	
	friend Point operator+(const Point & lhs,const Point &rhs);
	friend Point operator-(const Point & lhs,const Point &rhs);
	
private:
	//	坐标值
	double x;
	double y;
	double z;	//
};

#endif
