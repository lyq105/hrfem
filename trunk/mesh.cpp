


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
*	��ȡ��Ԫ���ͽڵ���
*
*/
	getline(element,buffer,'\n');
	nElem=atoi(buffer.c_str());

	getline(nodes,buffer,'\n');
	nNode=atoi(buffer.c_str());

//	�����ڴ�

	pNodes=new Nodes[nNode];
	pElements=new Elements[nElem];

/**
*
*	����ڵ���Ϣ 
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
*	����Ԫ��Ϣ
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
	��ʼ����һ��߽������ı�
  */
bool mesh2D::InitialBoundry1st(char* filename)
{
	return true;
};

/*!
	��ʼ���ڶ���߽������ı�
  */

bool mesh2D::InitialBoundry2nd(char* filename)
{
	return true;
};

/**
*
*	�ͷ�������ڴ�
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

