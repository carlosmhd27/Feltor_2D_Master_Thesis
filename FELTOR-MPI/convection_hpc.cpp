#include <iostream>
#include <iomanip>
#include <vector>

#ifdef TOEFL_MPI
#include <mpi.h>
#endif //FELTOR_MPI

#ifdef TOEFL_MPI
using HVec = dg::MHVec;
using DVec = dg::MDVec;
using DMatrix = dg::MDMatrix;
using IDMatrix = dg::MIDMatrix;
using Grid2d = dg::MPIGrid2d;
using CartesianGrid2d = dg::CartesianMPIGrid2d;
#define MPI_OUT if(rank==0)
#else //TOEFL_MPI
using HVec = dg::HVec;
using DVec = dg::DVec;
using DMatrix = dg::DMatrix;
using IDMatrix = dg::IDMatrix;
using Grid2d = dg::Grid2d;
using CartesianGrid2d = dg::CartesianGrid2d;
#define MPI_OUT
#endif //TOEFL_MPI

#include "dg/file/nc_utilities.h"

#include "convection.h"
#include "dg/algorithm.h"
#include "parameters.h"

int main( int argc, char* argv[])
{
  #ifdef TOEFL_MPI
      ////////////////////////////////setup MPI///////////////////////////////
      int provided;
      MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &provided);
      if( provided != MPI_THREAD_FUNNELED)
      {
          std::cerr << "wrong mpi-thread environment provided!\n";
          return -1;
      }
      int periods[2] = {false, true}; //non-, periodic
      int rank, size;
      MPI_Comm_rank( MPI_COMM_WORLD, &rank);
      MPI_Comm_size( MPI_COMM_WORLD, &size);

      int np[2];
      if(rank==0)
      {
          std::cin>> np[0] >> np[1];
          std::cout << "Computing with "<<np[0]<<" x "<<np[1]<<" = "<<size<<std::endl;
          assert( size == np[0]*np[1]);
      }
      MPI_Bcast( np, 2, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Comm comm;
      MPI_Cart_create( MPI_COMM_WORLD, 2, np, periods, true, &comm);
  #endif//TOEFL_MPI-

    //Parameter initialisation
    Json::Value js;
    Json::CharReaderBuilder parser;
    parser["collectComments"] = false;
    std::string errs;
    if( argc != 3)
    {
        std::cerr << "ERROR: Wrong number of arguments!\nUsage: "<< argv[0]<<" [inputfile] [outputfile]\n";
        return -1;
    }
    else
    {
        std::ifstream is(argv[1]);
        parseFromStream( parser, is, &js, &errs); //read input without comments
    }
    std::cout << js<<std::endl;
    const convection::Parameters p( js);
    p.display( std::cout);


    ////////////////////////////////set up computations///////////////////////////
    Grid2d grid( 0, p.lx, 0, p.ly, p.n, p.Nx, p.Ny, p.bc_x, p.bc_y);
    Grid2d grid_out( 0, p.lx, 0, p.ly, p.n_out, p.Nx_out, p.Ny_out, p.bc_x, p.bc_y);
    //create RHS
  	//exp is the explicit part and imp the implicit. Both defined as
  	//classes form the convection.h file (obviously)
    convection::ExplicitPart< CartesianGrid2d, DMatrix, DVec > exp( grid, p);
    convection::ImplicitPart< CartesianGrid2d, DMatrix, DVec > imp( grid, p.nu);
    //////////////////create initial vector///////////////////////////////////////
    dg::Gaussian g( p.posX*p.lx, p.posY*p.ly, p.sigma, p.sigma, p.amp); //gaussian width is in absolute values
    std::array<DVec,2> y0{
        dg::evaluate( g, grid), //n == Gaussian
        dg::evaluate(dg::zero,grid) // omega == 0
    };
    //////////////////initialisation of timekarniadakis and first step///////////////////
    double time = 0;
    dg::Karniadakis< std::array<DVec,2> > karniadakis( y0, y0[0].size(), p.eps_time);
    karniadakis.init( exp, imp, time, y0, p.dt);
    /////////////////////////////set up netcdf/////////////////////////////////////
    dg::file::NC_Error_Handle err;
    int ncid;
    err = nc_create( argv[2],NC_NETCDF4|NC_CLOBBER, &ncid);
    std::string input = js.toStyledString();
    err = nc_put_att_text( ncid, NC_GLOBAL, "inputfile", input.size(), input.data());
    int dim_ids[3], tvarID;

    err = dg::file::define_dimensions( ncid, dim_ids, &tvarID, grid_out);

    //Time ids
    int EtimeID, EtimevarID;
    err = dg::file::define_time( ncid, "energy_time", &EtimeID, &EtimevarID);

    //energy IDs
    std::array<int,4> invariantID, dissID;
    std::array<std::string,4> invariant_names { { "mass", "entropy", "kinetic", "curvature"}};
    std::array<std::string,4> dissipation_names{ { "mass_diss", "entropy_diss", "kinetic_diss", "curvature_diss"}};
    for( int i=0; i<4; i++)
    {
        err = nc_def_var( ncid, invariant_names[i].c_str(),  NC_DOUBLE, 1, &EtimeID, &invariantID[i]);
        err = nc_def_var( ncid, dissipation_names[i].c_str(), NC_DOUBLE, 1, &EtimeID, &dissID[i]);
    }
    err = nc_enddef(ncid);

    size_t Ecount[] = {1};
    size_t Estart[] = {0};

    err = nc_put_vara_double( ncid, EtimevarID, Estart, Ecount, &time);

    for( int i=0; i<4; i++)
    {
        err = nc_put_vara_double( ncid, invariantID[i], Estart, Ecount, &exp.invariants()[i]);
        err = nc_put_vara_double( ncid, dissID[i], Estart, Ecount, &exp.invariants_diffusion()[i]);
    }

    DVec transfer (dg::evaluate(dg::zero, grid));

    size_t count[3] = {1, grid_out.n()*grid_out.Ny(), grid_out.n()*grid_out.Nx()};
    size_t start[3] = {0, 0, 0}; // grid_out.n()*grid_out.Nx()
    int dataIDs[4];
    std::string names[4] = {"electrons", "ions", "potential", "vorticity"};
    IDMatrix interpolate = dg::create::interpolation( grid_out, grid);
    std::vector<DVec> transferD(4, dg::evaluate(dg::zero, grid_out)); // grid_out
    HVec transferH(dg::evaluate(dg::zero, grid_out)); // grid_out

    //field IDs

    for( unsigned i=0; i<4; i++){
        err = nc_def_var( ncid, names[i].data(), NC_DOUBLE, 3, dim_ids, &dataIDs[i]);}

    dg::blas2::symv( interpolate, y0[0],           transferD[0]);
    dg::blas2::symv( interpolate, y0[0],           transferD[1]);
    dg::blas2::symv( interpolate, exp.potential(), transferD[2]);
    dg::blas2::symv( interpolate, y0[1],           transferD[3]);

    for( int k=0;k<4; k++)
    {
        dg::assign( transferD[k], transferH);
        err = nc_put_vara_double( ncid, dataIDs[k], start, count, transferH.data() );
    }
    err = nc_put_vara_double( ncid, tvarID, start, count, &time);
    std::cout << "Exiting the fields" << std::endl;

    ///////////////////////////////////////Timeloop/////////////////////////////////
    const double mass0 = exp.invariants()[0], mass_blob0 = mass0 - grid.lx()*grid.ly();
    double E0 = exp.invariants()[1] + exp.invariants()[2], E1 = 0, diff = 0;
    dg::Timer t;
    t.tic();
    try
    {
#ifdef DG_BENCHMARK
    unsigned step = 0;
#endif //DG_BENCHMARK

    for( unsigned i=1; i<=p.maxout; i++)
    {

#ifdef DG_BENCHMARK
        dg::Timer ti;
        ti.tic();
#endif//DG_BENCHMARK

        for( unsigned j=0; j<p.itstp; j++)
        {
            karniadakis.step( exp, imp, time, y0);
            //store accuracy details
            {
                std::cout << "(m_tot-m_0)/m_0: "<< (exp.invariants()[0]-mass0)/mass_blob0<<"\t";
                E0 = E1;
                E1 = exp.invariants()[1] + exp.invariants()[2];
                diff = (E1 - E0)/p.dt;
                double diss = exp.invariants_diffusion()[1]+exp.invariants_diffusion()[2];
                std::cout << "diff: "<< diff<<" diss: "<<diss<<"\t";
                std::cout << "Accuracy: "<< 2.*(diff-diss)/(diff+diss)<<"\n";
            }

            Estart[0] += 1;
            // start[0] +=1;
            {
                std::cout << 0 <<std::endl;
                err = nc_open(argv[2], NC_WRITE, &ncid);

                err = nc_put_vara_double( ncid, EtimevarID, Estart, Ecount, &time);
                for( int i=0; i<4; i++)
                {
                    err = nc_put_vara_double( ncid, invariantID[i], Estart, Ecount, &exp.invariants()[i]);
                    err = nc_put_vara_double( ncid, dissID[i], Estart, Ecount, &exp.invariants_diffusion()[i]);
                }

                err = nc_close(ncid);
            }
        }
        //////////////////////////write fields////////////////////////
        start[0] = i;
        dg::blas2::symv( interpolate, y0[0], transferD[0]);
        dg::blas2::symv( interpolate, y0[0], transferD[1]);
        dg::blas2::symv( interpolate, exp.potential(), transferD[2]);
        dg::blas2::symv( interpolate, y0[1], transferD[3]);
        err = nc_open(argv[2], NC_WRITE, &ncid);
        for( int k=0;k<4; k++)
        {
            dg::assign( transferD[k], transferH);
            err = nc_put_vara_double( ncid, dataIDs[k], start, count, transferH.data() );
        }
        err = nc_put_vara_double( ncid, tvarID, start, count, &time);
        err = nc_close(ncid);

#ifdef DG_BENCHMARK
        ti.toc();
        step+=p.itstp;
        std::cout << "\n\t Step "<<step <<" of "<<p.itstp*p.maxout <<" at time "<<time;
        std::cout << "\n\t Average time for one step: "<<ti.diff()/(double)p.itstp<<"s\n\n"<<std::flush;
#endif//DG_BENCHMARK
    }
    }
    catch( dg::Fail& fail) {
        std::cerr << "CG failed to converge to "<<fail.epsilon()<<"\n";
        std::cerr << "Does Simulation respect CFL condition?\n";
    }
    t.toc();
    unsigned hour = (unsigned)floor(t.diff()/3600);
    unsigned minute = (unsigned)floor( (t.diff() - hour*3600)/60);
    double second = t.diff() - hour*3600 - minute*60;
    std::cout << std::fixed << std::setprecision(2) <<std::setfill('0');
    std::cout <<"Computation Time \t"<<hour<<":"<<std::setw(2)<<minute<<":"<<second<<"\n";
    std::cout <<"which is         \t"<<t.diff()/p.itstp/p.maxout<<"s/step\n";

    return 0;
}
