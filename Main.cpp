#include"IO.h"

int main()
{
	Schedule *sch = new Schedule();
	IO *io = new IO(sch);
	io->readInput();
	io->generateEdges();
	//system("pause");

	return 0;
}