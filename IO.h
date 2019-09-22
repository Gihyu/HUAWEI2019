#pragma once
#ifndef IO_H
#define IO_H
#include"Schedule.h"

typedef Eigen::Vector3d Point;

class IO
{
private:
	static Schedule * _schedule;


public:
	IO(Schedule * schedule);
	~IO();

	void readInput();
	void readSet();
	string GetStringFromCSV(string line, int nIdx);

	void generateEdges();

	double getDegAngle(Node* nd1, Node* nd2, Node* nd3);

};


#endif // !IO_H