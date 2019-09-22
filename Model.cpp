#include"Model.h"

Model::Model(Schedule* sch)
{
	_schedule = sch;
	_allNodeList = sch->getNodeList();
}


Model::~Model()
{
}

void Model::init()
{
	_model = IloModel(_env);
	_obj = IloAdd(_model, IloMinimize(_env));
	_solver = IloCplex(_model);


	vector<Node*> nodeList = _schedule->getNodeList();
	vector<Edge*> edgeList = _schedule->getAllEdgeList();

	int varINDEX = 0;
	_X_ij = IloNumVarArray(_env);
	_E_h_ij = IloNumVarArray(_env);
	_E_v_ij = IloNumVarArray(_env);
	for (int i = 0; i < nodeList.size(); i++)
	{
		vector<Node*> canVisitSet = nodeList[i]->getCanVisitSetByDis();
		for (int j = 0; j < canVisitSet.size(); j++)
		{
			string XvarName = "X_ij_" + to_string(nodeList[i]->getID()) + "_" + to_string(canVisitSet[j]->getID());
			_X_ij.add(IloNumVar(_env, 0, 1, ILOINT, XvarName.c_str()));

			string E_H_varName = "E_H_ij_" + to_string(nodeList[i]->getID()) + "_" + to_string(canVisitSet[j]->getID());
			_E_h_ij.add(IloNumVar(_env, 0, Util::Theta, ILOFLOAT, E_H_varName.c_str()));

			string E_V_varName = "E_V_ij_" + to_string(nodeList[i]->getID()) + "_" + to_string(canVisitSet[j]->getID());
			_E_v_ij.add(IloNumVar(_env, 0, Util::Theta, ILOFLOAT, E_V_varName.c_str()));

			_VarIndex[nodeList[i]->getID()][canVisitSet[j]->getID()] = varINDEX;

			varINDEX++;
		}
	}
	_model.add(_X_ij);
	_model.add(_E_h_ij);
	_model.add(_E_v_ij);

	_flow_blc = IloRangeArray(_env);  //X var balance cons
	for (int i = 0; i < nodeList.size(); i++)
	{
		IloExpr expr(_env);
		if (nodeList[i]->getNodeType() != SINK && nodeList[i]->getNodeType() != SOURCE)
		{
			vector<Edge*> outEdges = nodeList[i]->getOutEdgeSet();
			for (int oe = 0; oe < outEdges.size(); oe++)
			{
				int xVarIndex = _VarIndex[nodeList[i]->getID()][outEdges[oe]->getTailNode()->getID()];
				expr += _X_ij[xVarIndex];
			}

			vector<Edge*> inEdges = nodeList[i]->getInEdgeSet();
			for (int ie = 0; ie < inEdges.size(); ie++)
			{
				int xVarIndex = _VarIndex[inEdges[ie]->getHeadNode()->getID()][nodeList[i]->getID()];
				expr -= _X_ij[xVarIndex];
			}

			string consName = "flowBalance_SOURCE_Node_" + to_string(nodeList[i]->getID());
			_flow_blc.add(IloRange(_env, 0, expr, 0, consName.c_str()));
		}
		else if (nodeList[i]->getNodeType() == SOURCE)
		{
			vector<Edge*> outEdges = nodeList[i]->getOutEdgeSet();
			for (int oe = 0; oe < outEdges.size(); oe++)
			{	
				int xVarIndex = _VarIndex[nodeList[i]->getID()][outEdges[oe]->getTailNode()->getID()];
				expr += _X_ij[xVarIndex];
			}
			string consName = "flowBalance_Node_" + to_string(nodeList[i]->getID());
			_flow_blc.add(IloRange(_env, 1, expr, 1, consName.c_str()));
		}
		else if (nodeList[i]->getNodeType() == SINK)
		{
			vector<Edge*> inEdges = nodeList[i]->getInEdgeSet();
			for (int ie = 0; ie < inEdges.size(); ie++)
			{
				int xVarIndex = _VarIndex[inEdges[ie]->getHeadNode()->getID()][nodeList[i]->getID()];
				expr += _X_ij[xVarIndex];  //注意这里是加号
			}

			string consName = "flowBalance_SINK_Node_" + to_string(nodeList[i]->getID());
			_flow_blc.add(IloRange(_env, 1, expr, 1, consName.c_str()));
		}
		else 
		{
			cout << " BUG exists in initializing model !!!" << endl;
		}
		
		expr.end();
	 }
	_model.add(_flow_blc);
	 
	 _theta_h_cons = IloRangeArray(_env);
	 for (int i = 0; i < _X_ij.getSize(); i++)
	 {
		 IloExpr expr(_env);
		 expr -= _E_h_ij[i];
		 expr += Util::Theta*_X_ij[i];

		 string consName = "H_bound_by_Xij_" + to_string(i);
		 _theta_h_cons.add(IloRange(_env, 0, expr, IloInfinity, consName.c_str()));
		 expr.end();
	 }
	 _model.add(_theta_h_cons);

	 _theta_v_cons = IloRangeArray(_env);
	 for (int i = 0; i < _X_ij.getSize(); i++)
	 {
		 IloExpr expr(_env);
		 expr -= _E_v_ij[i];
		 expr += Util::Theta*_X_ij[i];

		 string consName = "V_bound_by_Xij_" + to_string(i);
		 _theta_v_cons.add(IloRange(_env, 0, expr, IloInfinity, consName.c_str()));
		 expr.end();
	 }
	 _model.add(_theta_v_cons);

	 _crNode_h_cons = IloRangeArray(_env);
	 for (int i = 0; i < nodeList.size(); i++)
	 {
		 IloExpr expr(_env);
		 if (nodeList[i]->getNodeType() != SINK && nodeList[i]->getNodeType() != SOURCE) //起终点没有这个约束
		 {
			 vector<Edge*> inEdges = nodeList[i]->getInEdgeSet();
			 for (int ie = 0; ie < inEdges.size(); ie++)
			 {
				 int xVarIndex = _VarIndex[inEdges[ie]->getHeadNode()->getID()][nodeList[i]->getID()];
				 expr -= _E_h_ij[xVarIndex];  //减的是h
			 }
			 if (nodeList[i]->getCorrectType() == VERTICAL)
			 {
				 expr += Util::Alpha_2; 
			 }
			 else if (nodeList[i]->getCorrectType() == HORIZONTAL)
			 {
				 expr += Util::Beta_2;
			 }
			 else
			 {
				 cout << " BUG exists in initializing model !!!" << endl;
			 }
			 string consName = "H_requirement_Node_" + to_string(nodeList[i]->getID());
			 _crNode_h_cons.add(IloRange(_env, 0, expr, IloInfinity, consName.c_str()));
		 }
		 expr.end();
	 }
	 _model.add(_crNode_h_cons);

	 _crNode_v_cons = IloRangeArray(_env);
	 for (int i = 0; i < nodeList.size(); i++)
	 {
		 IloExpr expr(_env);
		 if (nodeList[i]->getNodeType() != SINK && nodeList[i]->getNodeType() != SOURCE) //起终点没有这个约束
		 {
			 vector<Edge*> inEdges = nodeList[i]->getInEdgeSet();
			 for (int ie = 0; ie < inEdges.size(); ie++)
			 {
				 int xVarIndex = _VarIndex[inEdges[ie]->getHeadNode()->getID()][nodeList[i]->getID()];
				 expr -= _E_v_ij[xVarIndex];  //减的是v
			 }
			 if (nodeList[i]->getCorrectType() == VERTICAL)
			 {
				 expr += Util::Alpha_1;
			 }
			 else if (nodeList[i]->getCorrectType() == HORIZONTAL)
			 {
				 expr += Util::Beta_1;
			 }
			 else
			 {
				 cout << " BUG exists in initializing model !!!" << endl;
			 }
			 string consName = "V_requirement_Node_" + to_string(nodeList[i]->getID());
			 _crNode_v_cons.add(IloRange(_env, 0, expr, IloInfinity, consName.c_str()));
		 }
		 expr.end();
	 }
	 _model.add(_crNode_v_cons);


	 _correct_h_cons = IloRangeArray(_env);
	 for (int i = 0; i < nodeList.size(); i++)
	 {
		 vector<Node*> canVisitSet = nodeList[i]->getCanVisitSetByDis();
		 for (int j = 0; j < canVisitSet.size(); j++)
		 {
			 IloExpr expr(_env);

			 int varID = _VarIndex[nodeList[i]->getID()][canVisitSet[j]->getID()];
			 expr += _E_h_ij[varID];

			 if (nodeList[i]->getCorrectType() == VERTICAL)
			 {
				 vector<Edge*> inEdges = nodeList[i]->getInEdgeSet();  //compute E_h_ki
				 for (int k = 0; k < inEdges.size(); k++)
				 {	
					 //if (inEdges[k]->getHeadNode()->getID() != canVisitSet[j]->getID()) // BUG !!!!
					 //{
						 int ki_ID = _VarIndex[inEdges[k]->getHeadNode()->getID()][nodeList[i]->getID()];
						 expr -= _E_h_ij[ki_ID];
					 //}
					
				 }
			 }
			 expr -= _X_ij[varID]*Util::delta*nodeList[i]->getLinearDisFrom(canVisitSet[j]);  //Cij  delta  Xij
			 expr += 2 * Util::Theta*(1 - _X_ij[varID]);  // +2M(1-x)
			 string consName = "H_computation_From_Node_" + to_string(nodeList[i]->getID())+"_to_Node_"+to_string(canVisitSet[j]->getID());
			 _correct_h_cons.add(IloRange(_env, 0, expr, IloInfinity, consName.c_str()));
			 expr.end();
		 }
	 }
	 _model.add(_correct_h_cons);


	 _correct_v_cons = IloRangeArray(_env);
	 for (int i = 0; i < nodeList.size(); i++)
	 {
		 vector<Node*> canVisitSet = nodeList[i]->getCanVisitSetByDis();
		 for (int j = 0; j < canVisitSet.size(); j++)
		 {
			 IloExpr expr(_env);

			 int varID = _VarIndex[nodeList[i]->getID()][canVisitSet[j]->getID()];
			 expr += _E_v_ij[varID];

			 if (nodeList[i]->getCorrectType() == HORIZONTAL)
			 {
				 vector<Edge*> inEdges = nodeList[i]->getInEdgeSet();  //compute E_v_ki
				 for (int k = 0; k < inEdges.size(); k++)
				 {	
					 //if (inEdges[k]->getHeadNode()->getID() != canVisitSet[j]->getID())// BUG !!!!
					 //{
						 int ki_ID = _VarIndex[inEdges[k]->getHeadNode()->getID()][nodeList[i]->getID()];
						 expr -= _E_v_ij[ki_ID];
					 //}
				 }
			 }
			 expr -= _X_ij[varID] * Util::delta*nodeList[i]->getLinearDisFrom(canVisitSet[j]);  //Cij  delta  Xij
			 expr += 2 * Util::Theta*(1 - _X_ij[varID]);  // +2M(1-x)

			 string consName = "V_computation_From_Node_" + to_string(nodeList[i]->getID()) + "_to_Node_" + to_string(canVisitSet[j]->getID());
			 _correct_v_cons.add(IloRange(_env, 0, expr, IloInfinity, consName.c_str()));
			 expr.end();
		 }
	 }
	 _model.add(_correct_v_cons);

	 //OBJ!
	 IloExpr expr_obj(_env);
	 for (int i = 0; i < nodeList.size(); i++)
	 {
		 vector<Node*> canVisitSet = nodeList[i]->getCanVisitSetByDis();
		 for (int j = 0; j < canVisitSet.size(); j++)
		 {	
			 int x_varID = _VarIndex[nodeList[i]->getID()][canVisitSet[j]->getID()];
			 expr_obj += Util::w1*nodeList[i]->getLinearDisFrom(canVisitSet[j])*_X_ij[x_varID];
			 expr_obj += Util::w2*_X_ij[x_varID];
		 }
		
	 }
	 expr_obj -= Util::w2;
	 _obj.setExpr(expr_obj);
	 expr_obj.end();

}

void Model::solveMIP_arcModel()
{
	cout << endl << "---------------------------start to solve arcModel MIP-------------------------------" << endl;
	double t1 = clock();
	_solver = IloCplex(_model);
	string FileName = Util::OUTPUTPATH + string("arcModel") + string(".lp");
	_solver.exportModel(FileName.c_str());
	try
	{
		//_solver.setParam(IloCplex::EpGap, 0.002);
		_solver.setParam(IloCplex::TiLim, 600);
		//_solver.setParam(IloCplex::Param::Benders::Strategy,IloCplex::BendersFull);
		//_solver.setParam(IloCplex::RootAlg, IloCplex::Barrier);
		_solver.solve();

	}
	catch (IloException & e) { cerr << "Concert exception caught: " << e << endl; }
	catch (...) { cerr << "Unknown exception caught" << endl; }
	double t2 = clock();


	cout << endl << "Now print var size:" << endl;
	cout << "# Integer flow var X_ij:" << _X_ij.getSize() << endl;
	cout << "# H flow var E_H_ij:" << _E_h_ij.getSize() << endl;
	cout << "# V flow var E_V_ij:" << _E_v_ij.getSize() << endl;

	cout << "Now print cons size:" << endl;
	cout << "# x flow balance constraint:" << _flow_blc.getSize() << endl;
	cout << "# e_h bound by x_ij constraint:" << _theta_h_cons.getSize() << endl;
	cout << "# e_v bound by x_ij constraint:" << _theta_v_cons.getSize() << endl;
	cout << "# correction requirement for h computation:" << _crNode_h_cons.getSize() << endl;
	cout << "# correction requirement for v computation:" << _crNode_v_cons.getSize() << endl;
	cout << "# e_h calculation for ij:" << _correct_h_cons.getSize() << endl;
	cout << "# e_v calculation for ij:" << _correct_v_cons.getSize() << endl;

	cout << "Solution status: " << _solver.getStatus() << endl;

	if (_solver.getStatus() == IloAlgorithm::Infeasible)
	{
		cout << "!!!!!!!!!!!!!!!!!!!!!!!! Infeasibility in MIP sloving!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		return;
	}

	cout << "Optimal value: " << _solver.getObjValue() << endl;
	cout << "cplex running time is  " << t2 - t1 << "  ms." << endl << endl;

	//show details

	cout << "Now show var details" << endl << endl;
	IloNumArray xSoln = IloNumArray(_env);
	_solver.getValues(xSoln, _X_ij);

	IloNumArray hSoln = IloNumArray(_env);
	_solver.getValues(hSoln, _E_h_ij);

	IloNumArray vSoln = IloNumArray(_env);
	_solver.getValues(vSoln, _E_v_ij);
	
	cout << "x flow vars :" << endl;
	for (int i = 0; i < xSoln.getSize(); i++)
	{
		if (xSoln[i] > 0.001)
		{	
			int mIndex = 0;
			int nIndex = 0;
			for (int m = 0; m < _schedule->getNodeList().size(); m++)
			{	
				bool found = false;
				for (int n = 0; n < _schedule->getNodeList().size(); n++)
				{
					if (_VarIndex[m][n] == i)
					{
						mIndex = m;
						nIndex = n;
						found = true;
						break;
					}
				}
				if (found)
				{
					break;
				}
			}
			cout << _X_ij[i] << " : " << xSoln[i] << " and linearDis from node " << mIndex << " to " << nIndex << " is ";
			cout<< _allNodeList[mIndex]->getLinearDisFrom(_allNodeList[nIndex])<<endl;
		}
	}
	cout <<endl<< "h flow vars :" << endl;
	for (int i = 0; i < hSoln.getSize(); i++)
	{
		if (hSoln[i] > 0.001)
		{
			cout << _E_h_ij[i] << " : " << hSoln[i] << endl;
		}
	}
	cout << endl << "v flow vars :" << endl;
	for (int i = 0; i < vSoln.getSize(); i++)
	{
		if (vSoln[i] > 0.001)
		{
			cout << _E_v_ij[i] << " : " << vSoln[i] << endl;
		}
	}

}

void Model::end_model()
{
	_env.end();
}
