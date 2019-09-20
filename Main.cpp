#include"Model.h"
//#include"IO.h"

int main()
{
	Schedule *sch = new Schedule();
	IO *io = new IO(sch);
	io->readInput();
	io->generateEdges();
	//system("pause");

	Model* model = new Model(sch);
	model->init();
	model->solveMIP_arcModel();

	return 0;
}