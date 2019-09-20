#include "Node.h"

class OperNode
{
private:
	int _id;
	Node * _Node;

	OperNode * _frontOperNode;
	vector<OperNode *> _nextOperNodeSet;

	double _acmlDistance;  //acml: accumulated
	int _acmlCrNum;   //accumulated number of correction nodes

public:
	OperNode(int id, Node * nd) { _id = id, _Node = nd; }

	void setFrontOperNode(OperNode * opn) { _frontOperNode = opn; }
	OperNode * getFrontOperNode() { return _frontOperNode; }

	void pushNextOperNode(OperNode * opn) { _nextOperNodeSet.push_back(opn); }
	vector<OperNode *> getNextOPNset() { return _nextOperNodeSet; }
	void clearNextOPNset();

	void setAcmlDis(double acmldis){ _acmlDistance = acmldis;}
	double getAcmlDis() { return _acmlDistance; }

	void setAcmlCRnum(int acmlnum) { _acmlCrNum = acmlnum; }
	int getAcmlCRnum() { return _acmlCrNum; }

};
