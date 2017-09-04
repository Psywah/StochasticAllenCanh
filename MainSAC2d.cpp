#include "SacSolver.h"
using namespace dolfin;

/**
 * @brief main function
 *
 * @param argc
 * @param argv[]
 *
 * @return 
 */
int main()
{
    // load parameters
    Parameters para("user_defined_parameters");
    File para_file("./parameters.xml");
    para_file >> para;
    info(para, true);

    // generate mesh
    int Nc =(int)para["NumberOfCellPerDirection"];
    UnitSquareMesh mesh(Nc, Nc);
    
    // define variational solver
    info("initializing....");
    SacSolver sac(mesh, para);

    // solve
    info("solving....");
    sac.solve();


    return 0;
}
