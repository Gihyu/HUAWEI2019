#pragma once

#include"Node.h"
class Edge
{
private:
	int _id;
	int _avaiCount;
	Node * _headNode;
	Node * _tailNode;

	double _linearDis;

public:
	Edge(int id, Node*hn, Node*tn)
	{
		_id = id, _headNode = hn, _tailNode = tn, _linearDis = hn->getLinearDisFrom(tn);
	}

	void printEdgeInfo();

	double getLinearDis() { return _linearDis; }

	void setAvaiCount(int ac) { _avaiCount = ac; }
	int getAvaiCount() { return _avaiCount; }

	Node*getHeadNode(){return _headNode;}
	Node*getTailNode() { return _tailNode; }
};
