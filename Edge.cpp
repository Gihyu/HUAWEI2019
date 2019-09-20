#include"Edge.h"


void Edge::printEdgeInfo()
{
	cout << "Arc: ID " << _id << ", headNode is " << _headNode->getID() << ", tailNode is " << _tailNode->getID() << ", linearDistance is " << _linearDis << endl;
}