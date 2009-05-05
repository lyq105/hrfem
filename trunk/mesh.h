#ifndef _MESH_H_
#define _MESH_H_
#include "Point.h"


struct element2D{
	int i,j,k,l;
	int nNodes;
	Point nodes[4];
};
struct boundry1st2D
{
	int nNode;
	int Value;
};

struct boundry2nd2D
{
	int i,j;
	int nElem;
};

typedef Point Nodes;
typedef struct element2D Elements;
typedef struct boundry2nd2D Boundry2nds;
typedef struct boundry1st2D Boundry1sts;


class mesh2D
{
public:

	int Dimention;
	int nElem,nNode;
	
	
	Nodes* pNodes;
	Elements* pElements;
	Boundry1sts* pBoundry1st;
	Boundry2nds* pBoundry2nd;

	
	bool InitialMesh(char*,char*);
	bool InitialBoundry1st(char*);
	bool InitialBoundry2nd(char*);
	
	bool ReleaseMesh();
};



typedef mesh2D Mesh;






#endif
