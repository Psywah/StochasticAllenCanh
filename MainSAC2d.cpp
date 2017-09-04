#include "SacSolver.h"

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
    Parameters para("user_defined_parameters");
    File para_file("./parameters.xml");
    para_file >> para;
    info(para, true);

    int Nc =(int)para(["NumberOfCellPerDirection"]);
    UniteSquareMesh mesh(Nc, Nc);
    
    info("initializing....");
    SacSolver sac(mesh, para);
    info("solving....");
    sac.solve();


    return 0;
}
