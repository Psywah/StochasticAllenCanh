#ifndef __SAC_SOLVER_H
#define __SAC_SOLVER_H
#include <dolfin.h>


namespace dolfin
{
  /// This is a class for variational problems 
  class SacSolver
  {
  public:
    /// Constructor
    SacSolver(Mesh& mesh, Parameters& para);

    /// Destructor
    virtual ~SacSolver() {}

    void save_solution(Function& u);
    void save_solution(int sp, Function& u);
    double energy(Function& u);
    double spectrum(Function& u);
    void solve();
    //void load_solution(std::string str="backup_solution.xml");
    //void plot_solution();

    std::shared_ptr<Mesh> _mesh;

    // function space and forms
    std::shared_ptr<FunctionSpace> _V;
    std::shared_ptr<Form> _a;
    std::shared_ptr<Form> _L;
    PETScVector F;
    PETScMatrix A;

    // boundary conditions
    std::vector<std::shared_ptr<const DirichletBC>> bcs;

    // solutions
    std::shared_ptr<Function> _u;
    std::shared_ptr<Function> _u0;
    std::shared_ptr<Function> _uAverage;
    std::vector<std::shared_ptr<Function>> _us;
    std::vector<std::shared_ptr<Function>> _u0s;
    std::vector<std::shared_ptr<Function>> _uAverages;

    std::shared_ptr<File> _pvd_file_ave;
    std::shared_ptr<std::ofstream> _txt_file_ave;
    std::vector<std::shared_ptr<File>> _pvd_file;
    std::vector<std::shared_ptr<std::ofstream>> _txt_file;

    std::shared_ptr<Form> _energy;

#ifdef HAS_SLEPC
    // von compute spectrum
    std::shared_ptr<Form> _stiff;
    std::shared_ptr<Form> _mass;
    PETScMatrix Stiff, Mass;
    std::shared_ptr<SLEPcEigenSolver> _eigenSolver;
#endif

    Parameters para;
    std::shared_ptr<Constant> dw;
    double Save_dt, End_t, Dt;
    int repeat;

  };

}

#endif


