#include "IO.h"
Schedule * IO::_schedule = NULL;
IO::IO(Schedule * schedule)
{
	_schedule = schedule;
}

IO::~IO()
{
}

void IO::readInput()
{
	readSet();
}

void IO::readSet()
{
	string filename = Util::INPUTPATH + string("aftdata1.csv");
	ifstream file;
	file.open(filename.c_str());
	cout << "* Read location information from " << filename << "\n";
	string buff;
	char * token;
	char * tmp;
	getline(file, buff);
	getline(file, buff);
	while (getline(file, buff))
	{
		bool isVer = false;
		bool isHon = false;

		NODETYPE nodeType;
		CORRECTTYPE corctType;
		token = strtok_s((char *)buff.c_str(), ",", &tmp);
		int corID = atoi(token);
		token = strtok_s(NULL, ",", &tmp);
		double corX = atof(token);
		token = strtok_s(NULL, ",", &tmp);
		double corY = atof(token);
		token = strtok_s(NULL, ",", &tmp);
		double corZ = atof(token);
		token = strtok_s(NULL, ",", &tmp);
		//cout << string(token) << endl;
		if (string(token) == "A点")
		{
			nodeType = SOURCE;
			corctType = SOURCESINK;
		}
		else if (string(token) == "B点")
		{
			nodeType = SINK;
			corctType = SOURCESINK;
		}
		else
		{
			nodeType = NORMAL;
			int isVerInt = atoi(token);

			if (isVerInt)
			{
				corctType = VERTICAL;
			}
			else
			{
				corctType = HORIZONTAL;
			}
		}

		token = strtok_s(NULL, ",", &tmp);
		int isPro = atoi(token);

		Node * node = NULL;
		if (QUESTION_NUM == 1 || QUESTION_NUM == 2)
		{
			node = new Node(corID, corX, corY, corZ, 0, nodeType, corctType);
		}
		else
		{
			node = new Node(corID, corX, corY, corZ, isPro, nodeType, corctType);
		}
		_schedule->pushNode(node);
		//cout << corID << "," << corX << "," << corY<<","<<corZ << "," << isVer << "," << isHon << "," << isPro << endl;
	}

	cout << "----------------------------Print Node Informations---------------------------" << endl << endl;
	for (int i = 0; i < _schedule->getNodeList().size(); i++)
	{
		//cout<<_schedule->getNodeList()[i]->getID()<<endl;
		_schedule->getNodeList()[i]->printNodeInfo();
	}
	cout << "----------------------------Finish Node Informations---------------------------" << endl << endl;
	file.close();
}

string IO::GetStringFromCSV(string line, int nIdx)
{
	string str;
	int nSPos = 0;
	for (int i = 0; i < nIdx - 1; ++i)
	{
		nSPos = line.find(',', nSPos);
		++nSPos;
	}

	int nEPos = line.find(',', nSPos);
	if (nEPos != string::npos)
	{
		str = line.substr(nSPos, nEPos - nSPos);
	}
	else
	{
		str = line.substr(nSPos, line.size() - nSPos);
	}
	return str;
}


void IO::generateEdges()
{	
	vector<Node *> inputNodeSet = _schedule->getNodeList();
	int EdgID = 0;
	cout << "----------------------------Print All Edges Informations---------------------------" << endl << endl;
	for (int i=0; i < inputNodeSet.size()-1; i++)
	{
		for (int j = i + 1; j < inputNodeSet.size(); j++)
		{
			Edge * newEdg = new Edge(EdgID,inputNodeSet[i], inputNodeSet[j]);
			EdgID++;
			//newEdg->printEdgeInfo();
			_schedule->pushEdge_inAllset(newEdg);
			
			Edge * newEdg_2 = new Edge(EdgID, inputNodeSet[j], inputNodeSet[i]);  //反过来也要生成edge
			EdgID++;
			//newEdg_2->printEdgeInfo();
			_schedule->pushEdge_inAllset(newEdg_2);
		}
	}
	cout << "----------------------------Finish All Edges Informations---------------------------" << endl << endl;


	cout << "----------------------------Print Available Edges Informations---------------------------" << endl << endl;
	vector<Edge*> allEdges = _schedule->getAllEdgeList();

	int _avaiCount=0;
	for (auto & edg : allEdges)
	{	
		if (edg->getLinearDis()*Util::delta < Util::Theta)
		{	
			Node * headN = edg->getHeadNode();
			Node * tailN = edg->getTailNode();
			if (headN->getNodeType() != SINK && tailN->getNodeType()!=SOURCE)
			{
				_schedule->pushAvaiEdge(edg);
				edg->setAvaiCount(_avaiCount);
				_avaiCount++;

				headN->pushOutEdge(edg);
				headN->pushNodeInDisSet(tailN);  //前可以去后

				tailN->pushInEdge(edg);
				//edg->printEdgeInfo();
			}
		}
	}
	cout << " total available edges num " << _schedule->getAvaiEdges().size() << endl;
	cout << "----------------------------Finish Available Edges Informations---------------------------" << endl << endl;

	cout << "----------------------------Show Nodes can vist EdgeSetSize Informations---------------------------" << endl << endl;
	for (auto & nd : _schedule->getNodeList())
	{
		cout << "Node ID " << nd->getID() << ", canVistSize " << nd->getCanVisitSetByDis().size() << ", outEdgeSize " << nd->getOutEdgeSet().size() << ", inEdgeSize " << nd->getInEdgeSet().size() << endl;
	}
	cout << "----------------------------Finish Nodes can vist EdgeSetSize Informations---------------------------" << endl << endl;

	//check number of arcs and number of i can go to js
	int allcanGoTo = 0;
	for (auto * nod: _schedule->getNodeList())
	{
		allcanGoTo += nod->getCanVisitSetByDis().size();
	}
	cout << " total canVisitNodes num " << allcanGoTo << endl;


	/*vector<Edge*> outEdges = _schedule->getNodeList()[0]->getOutEdgeSet();
	for (int oe = 0; oe < outEdges.size(); oe++)
	{	
		int nextID = outEdges[oe]->getTailNode()->getID();
		cout << "hehe" << nextID <<endl;
	}*/
}
