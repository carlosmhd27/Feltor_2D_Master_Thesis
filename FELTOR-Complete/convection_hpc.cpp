#include <iostream>
#include <iomanip>
#include <vector>

#ifdef FELTOR_MPI
#include <mpi.h>
#endif //FELTOR_MPI

#include "dg/file/nc_utilities.h"
#include "dg/algorithm.h"
#include "convection.h"
#include "parameters.h"

#ifdef FELTOR_MPI
using HVec            = dg::MHVec;
using DVec            = dg::MDVec;
using DMatrix         = dg::MDMatrix;
using IDMatrix        = dg::MIDMatrix;
using Grid2d          = dg::MPIGrid2d;
using CartesianGrid2d = dg::CartesianMPIGrid2d;
#define MPI_OUT if(rank==0)
#else //FELTOR_MPI
using HVec            = dg::HVec;
using DVec            = dg::DVec;
using DMatrix         = dg::DMatrix;
using IDMatrix        = dg::IDMatrix;
using Grid2d          = dg::Grid2d;
using CartesianGrid2d = dg::CartesianGrid2d;
#define MPI_OUT
#endif //FELTOR_MPI


int main( int argc, char* argv[])
{
  #ifdef FELTOR_MPI
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
          np[0] = atoi(argv[3]);
          np[1] = atoi(argv[4]);
          std::cout << "Computing with "<<np[0]<<" x "<<np[1]<<" = "<<size<<std::endl;
          assert( size == np[0]*np[1]);
      }
      MPI_Bcast( np, 2, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Comm comm;
      MPI_Cart_create( MPI_COMM_WORLD, 2, np, periods, true, &comm);
      int num_args(5);
  #else
      int num_args(3);
  #endif//FELTOR_MPI-

    //Parameter initialisation
    Json::Value js;
    Json::CharReaderBuilder parser;
    parser["collectComments"] = false;
    std::string errs;

    if( argc != num_args)
    {
        #ifdef FELTOR_MPI
        MPI_OUT std::cerr << "ERROR: Wrong number of arguments!\nUsage: "<< argv[0]<<" [inputfile] [outputfile] [np0] [np1]\n";
        #else
        MPI_OUT std::cerr << "ERROR: Wrong number of arguments!\nUsage: "<< argv[0]<<" [inputfile] [outputfile]\n";
        #endif
        return -1;
    }
    else
    {
        std::ifstream is(argv[1]);
        parseFromStream( parser, is, &js, &errs); //read input without comments
    }
    MPI_OUT std::cout << js<<std::endl;
    const convection::Parameters p( js);
    MPI_OUT p.display( std::cout);


    ////////////////////////////////set up computations///////////////////////////
    Grid2d grid( 0, p.lx, 0, p.ly, p.n, p.Nx, p.Ny, p.bc_x_n, p.bc_y_n
                #ifdef FELTOR_MPI
                , comm
                #endif //FELTOR_MPI
    );
    Grid2d grid_out( 0, p.lx, 0, p.ly, p.n_out, p.Nx_out, p.Ny_out, p.bc_x_n, p.bc_y_n
                    #ifdef FELTOR_MPI
                    , comm
                    #endif //FELTOR_MPI
    );
    #ifndef FELTOR_MPI
    std::vector<unsigned> probes;
    std::vector<double>   probesx, probesy;
    if (p.save_pb){
        dg::Grid1d gridx( grid.x0(), grid.x1(), grid.n(), grid.Nx());
        dg::Grid1d gridy( grid.y0(), grid.y1(), grid.n(), grid.Ny());
        dg::HVec x_axis(dg::create::abscissas(gridx)), y_axis(dg::create::abscissas(gridy));
        for(unsigned i=0; i<p.probes.size(); i++){
            probes.push_back(p.probes[i][1] * p.Nx * p.n + p.probes[i][0]);
            probesx.push_back(x_axis[p.probes[i][0]]);
            probesy.push_back(y_axis[p.probes[i][1]]);
    }}
    #endif
    //create RHS
  	//exp is the explicit part and imp the implicit. Both defined as
  	//classes form the convection.h file (obviously)
    convection::ExplicitPart< CartesianGrid2d, DMatrix, DVec > exp( grid, p);
    //////////////////create initial vector///////////////////////////////////////
    dg::Gaussian g ( p.posX *p.lx, p.posY * p.ly, p.sigma,  p.sigma,  p.amp ); //gaussian width is in absolute values
    dg::TanhProfX tanh (p.x_b, p.tanh_width,  -1., p.profamp * (0.001 - p.nb), p.profamp * (p.nb)); // x_0, width, sign, B, A => B + A * 0.5(1 + sign tanh((x- x_0) / width))
    std::array<DVec,2> y0{
        dg::evaluate(tanh, grid),
        dg::evaluate(dg::zero,grid) // omega == 0
    };
    {DVec g2(dg::evaluate(g, grid));
    dg::blas1::axpby( 1., g2, 1., y0[0]);}
    //////////////////initialisation of the time stepper and first step///////////////////
    double time = 0, dt = p.dt;
    bool needs_to_init(true);
    dg::ExplicitMultistep< std::array<dg::DVec,2> > multi_stepper( "TVB-3-3", y0);
    dg::Adaptive<dg::ERKStep<std::array<dg::DVec, 2>>> adap_stepper( "Bogacki-Shampine-4-2-3", y0);
    if(p.Time_Step == "Multistep"){
        multi_stepper.init( exp, time, y0, p.dt);
        needs_to_init = false;
    }

    /////////////////////////////set up netcdf/////////////////////////////////////
    dg::file::NC_Error_Handle err, err_prb;
    int ncid, ncid_prb;
    MPI_OUT err     = nc_create( argv[2],NC_NETCDF4|NC_CLOBBER, &ncid);
    std::string nc_prb_fl = argv[2];
    nc_prb_fl = nc_prb_fl.insert(nc_prb_fl.size() - 3, "_prbs" );
    if(p.save_pb){
        MPI_OUT err_prb = nc_create( nc_prb_fl.c_str(),NC_NETCDF4|NC_CLOBBER, &ncid_prb);
    }
    std::string input = js.toStyledString();
    MPI_OUT err = nc_put_att_text( ncid, NC_GLOBAL, "inputfile", input.size(), input.data());
    int dim_ids[3], tvarID;

    MPI_OUT err = dg::file::define_dimensions( ncid, dim_ids, &tvarID, grid_out);

    // Probe IDs
    #ifndef FELTOR_MPI
    int probeID, EtimevarID, EtimeID;
    if (p.save_pb){
        //Time ids
        err_prb = dg::file::define_time( ncid_prb, "energy_time", &EtimeID, &EtimevarID);

        int probevarID, probevarxID, probevaryID;
        size_t count_probs[] = {probes.size()};
        size_t start_probs[] = {0};
        MPI_OUT err_prb = nc_def_dim( ncid_prb, "Probes",  probes.size(),  &probeID);
        MPI_OUT err_prb = nc_def_var( ncid_prb, "Probes",   NC_UINT,   1, &probeID, &probevarID);
        MPI_OUT err_prb = nc_put_vara_uint( ncid_prb, probevarID, start_probs, count_probs, probes.data());
        MPI_OUT err_prb = nc_def_var( ncid_prb, "Probes_x",  NC_DOUBLE, 1, &probeID, &probevarxID);
        MPI_OUT err_prb = nc_put_vara_double( ncid_prb, probevarxID, start_probs, count_probs, probesx.data());
        MPI_OUT err_prb = nc_def_var( ncid_prb, "Probes_y",  NC_DOUBLE, 1, &probeID, &probevaryID);
        MPI_OUT err_prb = nc_put_vara_double( ncid_prb, probevaryID, start_probs, count_probs, probesy.data());
    }
    #endif
    DVec transfer (dg::evaluate(dg::zero, grid));

    size_t start = 0, count = 1;
    int dataIDs[4];
    std::string names[4] = {"electrons", "ions", "potential", "vorticity"};
    IDMatrix interpolate = dg::create::interpolation( grid_out, grid);
    std::vector<DVec> transferD(4, dg::evaluate(dg::zero, grid_out)); // grid_out
    HVec transferH(dg::evaluate(dg::zero, grid_out)); // grid_out

    //field IDs

    for( unsigned i=0; i<4; i++){
        MPI_OUT err = nc_def_var( ncid, names[i].data(), NC_DOUBLE, 3, dim_ids, &dataIDs[i]);}

        dg::blas2::symv( interpolate, y0[0],           transferD[0]);
        dg::blas2::symv( interpolate, y0[0],           transferD[1]);
        dg::blas2::symv( interpolate, exp.potential(), transferD[2]);
        dg::blas2::symv( interpolate, y0[1],           transferD[3]);

    for( int k=0;k<4; k++)
    {
        dg::assign( transferD[k], transferH);
        dg::file::put_vara_double( ncid, dataIDs[k], start, grid_out, transferH);
    }
    MPI_OUT err = nc_put_vara_double( ncid, tvarID, &start, &count, &time);


    #ifndef FELTOR_MPI
    const unsigned prb_nmb = 4;
    std::array<std::vector<double>, prb_nmb> transfer_prb;
    int dim_prb_ids[2] = {EtimeID, probeID};
    size_t count_prb[2] = {1, probes.size()};
    size_t start_prb[2] = {0, 0};
    size_t Ecount[] = {1};
    size_t Estart[] = {0};
    int dataIDs_prb[prb_nmb];
    std::string names_prb[prb_nmb] = {"ions_probes", "potential_probes", "vorticity_probes", "vr_probes"};
    std::vector<dg::HVec> trnsfr_prbH(4, dg::evaluate(dg::zero, grid)); // grid_out

    if (p.save_pb){
        //Time
        MPI_OUT err_prb = nc_put_vara_double( ncid_prb, EtimevarID, Estart, Ecount, &time);
        dg::assign(y0[0],           trnsfr_prbH[0]);
        dg::assign(exp.potential(), trnsfr_prbH[1]);
        dg::assign(y0[1],           trnsfr_prbH[2]);
        dg::assign(exp.vradial(),   trnsfr_prbH[3]);
        //////////////first output ////////////
        for (unsigned k = 0; k < trnsfr_prbH.size(); k++){
            for (auto probe: probes){
                transfer_prb[k].push_back(trnsfr_prbH[k][probe]);
            }
            MPI_OUT err_prb = nc_def_var( ncid_prb, names_prb[k].data(),  NC_DOUBLE, 2,  dim_prb_ids, &dataIDs_prb[k]);
            MPI_OUT err_prb = nc_put_vara_double( ncid_prb, dataIDs_prb[k], start_prb, count_prb, transfer_prb[k].data());
        }
        MPI_OUT err_prb = nc_close(ncid_prb);
    }
    #endif
    ///////////////////////////////////////Timeloop/////////////////////////////////
    const double mass0 = exp.mass(), mass_blob0 = mass0 - grid.lx()*grid.ly();
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
            if(p.Time_Step == "Multistep" or (p.Time_Step == "Mixed" and time > p.tm_chng)){
                if(needs_to_init){
                    needs_to_init = false,
                    multi_stepper.init( exp, time, y0, p.dt);
                }
                multi_stepper.step( exp, time, y0);
                dt = p.dt;
            }
            // if(p.Time_Step == "Adaptive")
            else{
                adap_stepper.step( exp, time, y0, time, y0, dt, dg::pid_control, dg::l2norm, p.eps_time, p.eps_time);
            }
            //store accuracy details
            {
                MPI_OUT std::cout << "\t(m_tot-m_0)/m_0: "<< (exp.mass()-mass0)/mass_blob0<< " at time: " <<time <<"\n\t";
                MPI_OUT std::cout << "\tAt time: "<<time << " With dt: "<< dt <<"\n\t";
                MPI_OUT std::cout << "\n\t Step "<<step <<" of "<<p.itstp*p.maxout << "\n";
            }
            #ifndef FELTOR_MPI
            if (p.save_pb){
                start_prb[0] += 1;
                {
                MPI_OUT err_prb = nc_open(nc_prb_fl.c_str(), NC_WRITE, &ncid_prb);
                MPI_OUT err_prb = nc_put_vara_double( ncid_prb, EtimevarID, start_prb, Ecount, &time);
                dg::assign(y0[0],           trnsfr_prbH[0]);
                dg::assign(exp.potential(), trnsfr_prbH[1]);
                dg::assign(y0[1],           trnsfr_prbH[2]);
                dg::assign(exp.vradial(),   trnsfr_prbH[3]);
                for (unsigned k= 0; k < prb_nmb; k++){
                    for (unsigned l = 0; l < probes.size(); l++){
                        transfer_prb[k][l] = trnsfr_prbH[k][probes[l]];
                    }
                    MPI_OUT err_prb = nc_put_vara_double( ncid_prb, dataIDs_prb[k], start_prb, count_prb, transfer_prb[k].data());
                }
                MPI_OUT err_prb = nc_close(ncid_prb);}
            }
            #endif
        }
        //////////////////////////write fields/////////////////////////
        start = i;
        dg::blas2::symv( interpolate, y0[0], transferD[0]);
        dg::blas2::symv( interpolate, y0[0], transferD[1]);
        dg::blas2::symv( interpolate, exp.potential(), transferD[2]);
        dg::blas2::symv( interpolate, y0[1], transferD[3]);
        MPI_OUT err = nc_open(argv[2], NC_WRITE, &ncid);
        for( int k=0;k<4; k++)
        {
            dg::assign( transferD[k], transferH);
            dg::file::put_vara_double( ncid, dataIDs[k], start, grid_out, transferH);
        }
        MPI_OUT err = nc_put_vara_double( ncid, tvarID, &start, &count, &time);
        MPI_OUT err = nc_close(ncid);

#ifdef DG_BENCHMARK
        ti.toc();
        step+=p.itstp;
        MPI_OUT std::cout << "\n\t Step "<<step <<" of "<<p.itstp*p.maxout <<" at time "<<time;
        MPI_OUT std::cout << "\n\t Average time for one step: "<<ti.diff()/(double)p.itstp<<"s\n\n"<<std::flush;
#endif//DG_BENCHMARK
    }
    }
    catch( dg::Fail& fail) {
        MPI_OUT std::cerr << "CG failed to converge to "<<fail.epsilon()<<"\n";
        MPI_OUT std::cerr << "Does Simulation respect CFL condition?\n";
    }
    t.toc();
    unsigned hour = (unsigned)floor(t.diff()/3600);
    unsigned minute = (unsigned)floor( (t.diff() - hour*3600)/60);
    double second = t.diff() - hour*3600 - minute*60;
    MPI_OUT std::cout << std::fixed << std::setprecision(2) <<std::setfill('0');
    MPI_OUT std::cout <<"Computation Time \t"<<hour<<":"<<std::setw(2)<<minute<<":"<<second<<"\n";
    MPI_OUT std::cout <<"which is         \t"<<t.diff()/p.itstp/p.maxout<<"s/step\n";

    return 0;
}
