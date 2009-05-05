


#include "mesh.h"

#include <iostream>
#include <fstream>
#include <string>

using std::ifstream;
using std::string;
using std::getline;

bool mesh2D::InitialMesh(char* ElementFilename,char* NodesFileName)		
{

	ifstream element(ElementFilename);
	ifstream nodes(NodesFileName);

	string buffer;
/*!
*
*	读取单元数和节点数
*
*/
	getline(element,buffer,'\n');
	nElem=atoi(buffer.c_str());

	getline(nodes,buffer,'\n');
	nNode=atoi(buffer.c_str());

//	分配内存

	pNodes=new Nodes[nNode];
	pElements=new Elements[nElem];

/**
*
*	处理节点信息 
*
*/
	for (int i=0;i<nNode;i++)
	{
		getline(nodes,buffer,'\n');
		if (buffer.substr(0,1)=="#" )
		{
			getline(nodes,buffer,'\n');
		}
		if (buffer.empty())
		{
			getline(nodes,buffer,'\n');
		}
		pNodes[i].setx(atof(buffer.substr(8,20).c_str()));
		pNodes[i].sety(atof(buffer.substr(28,20).c_str()));
	}

/**
*
*	处理单元信息
*
*/	
	for (i=0;i<nElem;i++)
	{
		getline(element,buffer,'\n');
		if (buffer.substr(0,1)=="#")
		{
			getline(element,buffer,'\n');
		}
		if (buffer.empty())
		{
			getline(element,buffer,'\n');
		}
		pElements[i].i=atoi(buffer.substr(0,6).c_str());
		pElements[i].j=atoi(buffer.substr(6,6).c_str());
		pElements[i].k=atoi(buffer.substr(12,6).c_str());
		pElements[i].l=atoi(buffer.substr(18,6).c_str());
	}

	nodes.close();
	element.close();

	return true;
}


/*!
	初始化第一类边界条件的边
  */
bool mesh2D::InitialBoundry1st(char* filename)
{
	return true;
};

/*!
	初始化第二类边界条件的边
  */

bool mesh2D::InitialBoundry2nd(char* filename)
{
	return true;
};

/**
*
*	释放网格的内存
*
*/
bool mesh2D::ReleaseMesh()
{

	delete [] pElements;
	delete [] pNodes;
	delete [] pBoundry1st;
	delete [] pBoundry2nd;

	return true;
}

