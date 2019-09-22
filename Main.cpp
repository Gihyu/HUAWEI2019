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
	//model->init();
	model->init_Q3_nonPro();
	//model->init_forSet1_250_340();
	model->solveMIP_arcModel();

	//cout << io->getDegAngle(sch->getNodeList()[0], sch->getNodeList()[1], sch->getNodeList()[2]) << endl;


	return 0;
}