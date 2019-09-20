#pragma once
#include"Edge.h"
class Schedule
{
private:
	vector<Node*> _nodeList;
	vector<Edge*> _allEdgeList;

	vector<Edge*> _availableEdgeList;  //linearDistan < delta*theta.  they will be showed in class Node, inEdges or outEdges.

public:
	Schedule();
	~Schedule();

	vector<Node*> getNodeList() { return _nodeList; }
	void setNodeList(vector<Node*> nodeList) { _nodeList = nodeList; }
	void pushNode(Node * node) { _nodeList.push_back(node); }

	vector<Edge*> getAllEdgeList() { return _allEdgeList; }
	void setAllEdgeList(vector<Edge*> edglist) { _allEdgeList = edglist; }
	void pushEdge_inAllset(Edge * edg) { _allEdgeList.push_back(edg); }

	vector<Edge*> getAvaiEdges() { return _availableEdgeList; }
	void setAvaiEdges(vector<Edge*> edglist) { _availableEdgeList = edglist; }
	void pushAvaiEdge(Edge * edg) { _availableEdgeList.push_back(edg); }
};