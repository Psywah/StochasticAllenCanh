#include "SacSolver.h"
#include "AllenCahn2D.h"
#include "EnergyForm.h"
#include "SpectrumForm.h"
#include <iostream>
#include <fstream>
#include <math.h>


using namespace dolfin;

// Initial conditions
class InitialConditions : public Expression
{
public:
  void eval(Array<double>& values, const Array<double>& x) const
  {
    double d = sqrt((x[0]-0.5)*(x[0]-0.5)+(x[1]-0.5)*(x[1]-0.5)) - 0.25;
   // values[0] = tanh(d/sqrt(2)/EPS);
   values[0] = tanh(d/(sqrt(2)*EPS));
  }
};

void SacSolver::save_solution(int sp, Function& u);
{
    *(_pvd_file[sp]) << u;
    for (int i=0; i< u.vector()->size(); ++i) 
    {
        *(_txt_file[sp])  << (*u->vector())[i] << std::endl;
    }
}

void SacSolver::save_solution(Function& u)
{
    *(_pvd_file_ave) << u;
    for (int i=0; i< u.vector()->size(); ++i) 
    {
        *(_txt_file_ave)  << (*u->vector())[i] << std::endl;
    }

}

double SacSolver::energy(Function& u)
{
    *(_u->vector()) = *(u.vector());
    return assemble(_energy);
}

double SacSolver::spectrum(Function& u)
{

#ifdef HAS_SLEPC
    *(_u->vector()) = *(u.vector());
    assemble(Stiff, *_stiff);
    _eigenSolver->parameters["spectrum"]= "smallest magnitude";
    _eigenSolver->solve();
    double r, c;
    PETScVector rx, cx;
    _eigenSolver->get_eigenpair(r, c, rx, cx);
    return r;
    //info("Smallest eigenvalue: %f ", r); 
#else
    info("please install dolfin with slepc");
    return 0.0;
#endif
}

void SacSolver::solve()
{
    double t = Dt;
    double nextSaveTime = t;

    while(t < End_t +DOLFIN_EPS)
    {
        begin("Step forwad");
        info("Time %.4f", t);
        begin("random repeat");
        _uAverage->vector()->zero();
        for(int i = 0; i <repeat; ++i)
        {
            info("Stochostic Process %d", i);

            *(_u0->vector()) = *(_u0s[i]->vector());
            dw = dolfin:rand() - .5;

            assemble(F,*_L);
            solve(*_A, *(_u->vector()), F);
            *(_u0s[i]->vector()) = *(_u-vector());
            *(_uAverage->vector()) += *(_u->vector());

            double e = energy(*_u);
            double lambda = specturm(*_u);
            
            info("energy: %f, specturm: %f", e, lambda);
            if(abs(t - nextSaveTime) <1.e-5)
                    save_solution(i, *_u);
        }
        *(_uAverage->vector()) /= repeat;
        double e = energy(*_uAverage);
        double lambda = specturm(*_uAverage);
        info("Average energy: %f, Average specturm: %f", e, lambda);
        if(abs(t - nextSaveTime) <1.e-5)
        {
            save_solution(*_uAverage);
            nextSaveTime = t + dt * floor(Save_dt /Dt);
            if(nextSaveTime > End_T)
                nextSaveTime = End_t;
        }
        end();
        end();

        t = t +Dt;
    }

    info("solved!");

}

SacSolver::SacSolver( Mesh& mesh, Parameters& _para): para(_para)
{
    double Eps = (double)para(["epsilon"]),
           Sigma = (double)para(["sigma"]);
    Save_dt = (double)para(["SaveDt"]);
    End_t = (double)para(["EndT"]);
    Dt = (double)para(["TimeStep"]);
    repeat = (int)para(["repeat"]);

    auto epsilon = std::make_shared<Constant>(Eps);
    auto dt = std::make_shared<Constant>(Dt); 
    auto sigma = std::make_shared<Constant>(Sigma);


    dolfin::seed(777);
    dw = std::make_shared<Constant>(dolfin::rand() - .5);
   
    _mesh = reference_to_no_delete_pointer(mesh);
    _V = std::make_shared<AllenCahn2D::FunctionSpace>(_mesh);
    int num_v = mesh.num_vertices(), num_c = mesh.num_cells();
    info("# of Vertices %d, # of Cells %d, # of Unknowns %d", dolfin::MPI::sum(MPI_COMM_WORLD,num_v), dolfin::MPI::sum(MPI_COMM_WORLD,num_c), _V->dim());


    _u = std::make_shared<Function>(_V);
    _u0 = std::make_shared<Function>(_V);
    _uAverage = std::make_shared<Function>(_V);
    InitialConditions u_initial;
    _u0->interpolate(u_initial);
    for(int i = 0; i<repeat; ++i)
    {
        std::shared_ptr<Function> tmp1 = std::make_shared<Function>(_V);
        std::shared_ptr<Function> tmp2 = std::make_shared<Function>(_V);
        std::shared_ptr<Function> tmp3 = std::make_shared<Function>(_V);
        *(tmp1->vector()) = *(_u0->vector());
        _u0s.push_back(tmp1);
        _us.push_back(tmp3);
        _uAverages.push_back(tmp2);

        std::shared_ptr<File> tmpfile = std::make_shared<File>(std::string("./result/sp") + std::to_string(repeat - i) + std::string("/allen_cahn.pvd"));
        _pvd_file.push_back(tmpfile);
        std::shared_ptr<std::ofstream> tmpfileTxt = std::make_shared<std::ofstream>();
        std::string filename = std::string("./result/sp") +std::to_string(repeat - i) + std::string("/phi.txt");
        tmpfileTxt->open(filename, std::iso_base::app);
        _txt_file.push_back(tmpfileTxt);
    }
    _pvd_file_ave = std::make_shared<File>(std::string("./result/ave") +  std::string("/allen_cahn.pvd"));
    _txt_file_ave = std::make_shared<std::ofstream>();
    _txt_file_ave->open("./result/ave/phi.txt", std::iso_base::app);

    std::map<std::string, std::shared_ptr<const GenericFunction>> coef_list
                = { {"u",       _u},
                    {"u0",      _u},
                    {"dw",      dw},
                    {"eps",     epsilon},
                    {"dt",      dt},
                    {"sigma",   sigma}};
    info("number of coefficients %d", coef_list.size());

    _a = std::make_shared<AllenCahn2D::bilinearForm>(_V,_V);
    _L = std::make_shared<AllenCahn2D::LinearForm>(_V);
    _energy = std::makr_shared<EnergyForm::Functional>(_mesh);

    _a->set_some_coefficients(coef_list);
    _L->set_some_coefficients(coef_list);
    _energy->set_some_coefficients(coef_list);

    assemble(A, *_a);
    assemble(F, *_L);

#ifdef HAS_SLEPC
    info("has slepc");
    _stiff =std::make_shared<SpectrumForm::Form_a>(_V,_V);
    _mass =std::make_shared<SpectrumForm::Form_m>(_V,_V);
    _stiff->set_some_coefficients(coef_list);
    _mass->set_some_coefficients(coef_list);
    assemble(Stiff, *_stiff);
    assemble(Mass, *_mass);
    _eigenSolver = std::make_shared<SLEPcEigenSolver>(reference_to_no_delete_pointer(Stiff), reference_to_no_delete_pointer(Mass));
    _eigenSolver->parameters["spectrum"]= "smallest magnitude";
#else
    info("please install dolfin with slepc");
#endif

}


