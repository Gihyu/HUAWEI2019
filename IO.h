#pragma once
#ifndef IO_H
#define IO_H
#include"Schedule.h"
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

};


#endif // !IO_H