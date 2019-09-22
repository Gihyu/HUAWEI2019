#pragma once
#include"IO.h"

class Model
{
private:
	Schedule* _schedule;
	vector<Node* >_allNodeList;
	
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

	IloRangeArray _Q3_p_cons;

	int _VarIndex[NODE_NUM][NODE_NUM];

	// Q3
	IloNumVarArray _Z_h_ij;
	IloNumVarArray _Y_h_ij;

	IloRangeArray H_Y_Z_BOUND_1;
	IloRangeArray H_Y_Z_BOUND_2;
	IloRangeArray H_Y_Z_BOUND_3;
	IloRangeArray H_Y_Z_BOUND_4;
	IloRangeArray H_Y_Z_BOUND_5;
	IloRangeArray H_Y_Z_BOUND_6;

	IloNumVarArray _Z_v_ij;
	IloNumVarArray _Y_v_ij;

	IloRangeArray V_Y_Z_BOUND_1;
	IloRangeArray V_Y_Z_BOUND_2;
	IloRangeArray V_Y_Z_BOUND_3;
	IloRangeArray V_Y_Z_BOUND_4;
	IloRangeArray V_Y_Z_BOUND_5;
	IloRangeArray V_Y_Z_BOUND_6;

public:
	Model(Schedule *sch);
	~Model();
	void init();
	void init_Q3_nonPro();
	void end_model();

	void init_forSet1_250_340();

	void solveMIP_arcModel();
};