#pragma once
#include"Schedule.h"

typedef	IloArray<IloArray<IloNumVarArray>> NumVar3Matrix;
typedef IloArray<IloArray<IloNumArray>> Num3Matrix;

typedef IloArray<IloNumVarArray> NumVar2Matrix;
typedef IloArray<IloNumArray> Num2Matrix;

class Model
{
private:
	NumVar2Matrix X_ij;
	NumVar2Matrix E_h_ij;
	NumVar2Matrix E_v_ij;

public:
	
};