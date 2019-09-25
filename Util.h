﻿#pragma once

#include<iostream>
using namespace std;
#include<stdio.h>
#include <algorithm>
#include<string>
#include<ctime>
#include<vector>
#include<stack>
#include<fstream>
#include<sstream>
#include<ctime>
#include<process.h>
#include<random>
#include<array>
#include<map>
#include <iomanip>
#include<ilcplex/ilocplex.h>
#include <Eigen/Dense>

//#define NODE_NUM 613
#define NODE_NUM 327
#define QUESTION_NUM 1
#define INPUTSET_NUM 2

class Util
{
public:
	static string INPUTPATH;
	static string OUTPUTPATH;

	static double delta;

	static double Theta;
	static double Alpha_1;
	static double Alpha_2;
	static double Beta_1;
	static double Beta_2;

	static double w1;
	static double w2;

};


