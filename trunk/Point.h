#ifndef _POINT_H_
#define _POINT_H_

class Point
{
public:
	Point(Point&);
	Point(double =0,double =0,double =0);
	~Point(){};
public:
	//	��������ֵ
	inline double getx() {return x;}
	inline double gety() {return y;}
	inline double getz() {return z;}
	//  ����const����
	inline double getx()const{return x;}
	inline double gety()const{return y;}
	inline double getz()const{return z;}
	
	//	ȷ������ֵ
	inline void setx(double a){ x = a;}
	inline void sety(double b){ y = b;}
	inline void setz(double c){ z = c;}
	//	���������
	Point &operator = (const Point & pt);
	Point &operator +=(const Point & pt);
	
	friend Point operator * (Point& pt,double rhs);
	friend Point operator * (double lhs,Point &pt);	
	friend Point operator+(const Point & lhs,const Point &rhs);
	friend Point operator-(const Point & lhs,const Point &rhs);
	
private:
	//	����ֵ
	double x;
	double y;
	double z;	//
};

#endif
