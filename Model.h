#pragma once
#include"IO.h"

class Model
{
private:
	Schedule* _schedule;
	
	IloModel _model;
	IloEnv _env;
	IloCplex _solver;
	IloObjective _obj;

	IloNumVarArray _X_ij;
	IloNumVarArray _E_h_ij;
	IloNumVarArray _E_v_ij;

	IloRangeArray _flow_blc;
	IloRangeArray _theta_h_cons;
	IloRangeArray _theta_v_cons;
	IloRangeArray _crNode_h_cons;
	IloRangeArray _crNode_v_cons;
	IloRangeArray _correct_h_cons;
	IloRangeArray _correct_v_cons;

	int _VarIndex[NODE_NUM][NODE_NUM];


public:
	Model(Schedule *sch);
	~Model();
	void init();
	void end_model();

	void solveMIP_arcModel();
};