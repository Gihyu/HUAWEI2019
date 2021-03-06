#include "Node.h"


bool Node::canVisitThisByDis(Node * nd)
{
	if (_canVisitSetByDis.empty())
	{
		return false;
	}
	else
	{
		for (auto& canVist : _canVisitSetByDis)
		{
			if (canVist == nd)
			{
				return true;
			}
		}
	}
	return false;
}

double Node::getLinearDisFrom(Node* nd)
{
	double nextX = nd->getX();
	double nextY = nd->getY();
	double nextZ = nd->getZ();

	//return sqrt(pow(_x - nextX, 2) + pow(_y - nextY, 2) + pow(_z - nextZ, 2));
	return sqrt((_x - nextX)*(_x - nextX) + (_y - nextY)*(_y - nextY) + (_z - nextZ)*(_z - nextZ));
}

void Node::printNodeInfo()
{
	cout << "Node: ID " << _id << ", X " << _x << ", Y " << _y << ", Z " << _z << ", Type " << _nodeType << ", CrType " << _crType << ", canVisitListSize " << _canVisitSetByDis.size() <<" isProblemNode "<< _mayBePerturbed<< endl;
}

bool Node::isProblemNode()
{
	if (_mayBePerturbed == 1)
	{
		return true;
	}
	else
	{
		return false;
	}
}