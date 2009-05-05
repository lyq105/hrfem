/**********************************************************
 *
 *	Geometry.h
 *	The interface of Geomtry.cpp
 *
 *	 Created: March 13,2009 22:11
 *	Modified: Marth 13,2009 22:11
 * *******************************************************/



#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_


// Data Structure of Domain's Geometry;


//node
typedef struct{
	double x,y,z;
	int dim=1;
} node;

//element 
typedef struct{
	int i,j,k,l,m,n,p,q; // Element's node number; 
	int mertrial;
} element;

//side
typedef struct{
	int i,j; 	// The ends of a side;
	int mark; 	// Is domain's boundry;
}side;
typedef struct{
	int nNodes;
	int nElements;
	int DimentionofDomain;
	node*  pNodes;
	element* pElement;
}

// Methods of Using Geometry Data.

bool ReadGeometryData(char* fname);


#endif //_GEOMETRY_H_
