#include  "SolveSequence.h"
#include  "Equation.h"


// this is the main function of this project

int main()
{

	FEM fe;

	fe.Initial("element","node");
	fe.FormMatrix();
	return 0;

}
