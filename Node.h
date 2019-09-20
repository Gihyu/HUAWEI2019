#pragma once

#include"Util.h"
class OperNode;
class Edge;

enum NODETYPE
{
	NORMAL = 0,
	SOURCE = 1,
	SINK = 2,
};

enum CORRECTTYPE
{
	HORIZONTAL = 0,
	VERTICAL = 1,
	SOURCESINK = 2,
};

class Node
{
private:
	int _id;

	double _x;   //coordinate
	double _y;
	double _z;

	NODETYPE _nodeType;
	CORRECTTYPE _crType;

	int _mayBePerturbed;   //1, this node may not success in correcting with question 2;   0,otherwise
						   // 1 means this node might fail with probility 20%

	int _situation;  // mark the situation of the Block (Dijkstra? BFS?) in order to run faster
					 // 0 means notFound ; 1 means Searching ; 2 means visted
					 // initial is 0

	vector<Node* > _canVisitSetByDis;   //initialized by the maximum linear distance, e.g., 30 in our competition

	vector<Edge* > _inEdges;
	vector<Edge* > _outEdges;

public:
	Node(int id, double x, double y, double z, int problem,NODETYPE nodeType,CORRECTTYPE crType)
	{
		_id = id, _x = x ,_y = y,_z = z, _mayBePerturbed = problem,_nodeType = nodeType,_crType = crType;
	}
	int getID() { return _id; }
	double getX() { return _x; }
	double getY() { return _y; }
	double getZ() { return _z; }

	void setSituation( int st ){_situation = st; }
	int getSituation() { return _situation; }

	void pushNodeInDisSet(Node * nd) { _canVisitSetByDis.push_back(nd); }
	vector<Node* > getCanVisitSetByDis(){ return _canVisitSetByDis; }
	bool canVisitThisByDis(Node * nd);

	double getLinearDisFrom(Node *nd);

	void setNodeType(NODETYPE ndtype) { _nodeType = ndtype; }
	NODETYPE getNodeType() { return _nodeType; }

	void setCorrectType(CORRECTTYPE crtype) { _crType = crtype; }
	CORRECTTYPE getCorrectType() { return _crType; }

	void printNodeInfo();

	void pushInEdge(Edge* edg) { _inEdges.push_back(edg); }
	vector<Edge*> getInEdgeSet() { return _inEdges; }

	void pushOutEdge(Edge* edg) { _outEdges.push_back(edg); }
	vector<Edge*> getOutEdgeSet() { return _outEdges; }
	
};