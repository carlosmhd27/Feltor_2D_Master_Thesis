#include <iostream>
#include <iomanip>
#include <vector>

#include "dg/file/nc_utilities.h"

#include "convection.h"
#include "dg/algorithm.h"
#include "parameters.h"

int main( int argc, char* argv[])
{
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
    dg::Grid2d grid( 0, p.lx, 0, p.ly, p.n, p.Nx, p.Ny, p.bc_x, p.bc_y);
    dg::Grid2d grid_out( 0, p.lx, 0, p.ly, p.n_out, p.Nx_out, p.Ny_out, p.bc_x, p.bc_y);
    std::vector<unsigned> probes;
    std::vector<double>   probesx, probesy;
    if (p.save_pb){
    dg::Grid1d gridx( grid.x0(), grid.x1(), grid.n(), grid.Nx());
    dg::Grid1d gridy( grid.y0(), grid.y1(), grid.n(), grid.Ny());
    dg::HVec x_axis(dg::create::abscissas(gridx)), y_axis(dg::create::abscissas(gridy));
    for(unsigned i=0; i<p.probes.size(); i++){
        probes.push_back(p.probes[i][1] * p.n * p.Ny * p.n + p.probes[i][0] * p.n);
        probesx.push_back(x_axis[p.n * p.probes[i][0]]);
        probesy.push_back(y_axis[p.n * p.probes[i][1]]);
        std::cout << "Xpositions" << probesx[i] << '\n';
    }}
    //create RHS
	//exp is the explicit part and imp the implicit. Both defined as
	//classes form the convection.h file (obviously)
    convection::ExplicitPart< dg::CartesianGrid2d, dg::DMatrix, dg::DVec > exp( grid, p);
    convection::ImplicitPart< dg::CartesianGrid2d, dg::DMatrix, dg::DVec > imp( grid, p.nu);
    //////////////////create initial vector///////////////////////////////////////
    dg::Gaussian g( p.posX*p.lx, p.posY*p.ly, p.sigma, p.sigma, p.amp); //gaussian width is in absolute values
    std::array<dg::DVec,2> y0{
        dg::evaluate( g, grid), //n == Gaussian
        dg::evaluate(dg::zero,grid) // omega == 0
    };
    //////////////////initialisation of timekarniadakis and first step///////////////////
    double time = 0;
    dg::Karniadakis< std::array<dg::DVec,2> > karniadakis( y0, y0[0].size(), p.eps_time);
    karniadakis.init( exp, imp, time, y0, p.dt);
    /////////////////////////////set up netcdf/////////////////////////////////////
    dg::file::NC_Error_Handle err;
    int ncid;
    err = nc_create( argv[2],NC_NETCDF4|NC_CLOBBER, &ncid);
    std::string input = js.toStyledString();
    err = nc_put_att_text( ncid, NC_GLOBAL, "inputfile", input.size(), input.data());

    //Fields dimenstions
    int dim_ids[3], tvarID;
    err = dg::file::define_dimensions( ncid, dim_ids, &tvarID, grid_out);

    //Time ids
    int EtimeID, EtimevarID;
    err = dg::file::define_time( ncid, "energy_time", &EtimeID, &EtimevarID);

    // Probe IDs
    int probeID;
    if (p.save_pb){
        int probevarID, probevarxID, probevaryID;
        size_t count_probs[] = {probes.size()};
        size_t start_probs[] = {0};
        err = nc_def_dim( ncid, "Probes",  probes.size(),  &probeID);
        err = nc_def_var( ncid, "Probes",   NC_UINT,   1, &probeID, &probevarID);
        err = nc_put_vara_uint( ncid, probevarID, start_probs, count_probs, probes.data());
        err = nc_def_var( ncid, "Probes_x",  NC_DOUBLE, 1, &probeID, &probevarxID);
        err = nc_put_vara_double( ncid, probevarxID, start_probs, count_probs, probesx.data());
        err = nc_def_var( ncid, "Probes_y",  NC_DOUBLE, 1, &probeID, &probevaryID);
        err = nc_put_vara_double( ncid, probevaryID, start_probs, count_probs, probesy.data());
    }
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

    dg::DVec transfer (dg::evaluate(dg::zero, grid));
    size_t count[3] = {1, grid_out.n()*grid_out.Ny(), grid_out.n()*grid_out.Nx()};
    size_t start[3] = {0, 0, 0}; // grid_out.n()*grid_out.Nx()
    int dataIDs[4];
    std::string names[4] = {"electrons", "ions", "potential", "vorticity"};
    dg::IDMatrix interpolate = dg::create::interpolation( grid_out, grid);
    std::vector<dg::DVec> transferD(4, dg::evaluate(dg::zero, grid_out)); // grid_out
    dg::HVec transferH(dg::evaluate(dg::zero, grid_out)); // grid_out
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


    // std::vector<dg::DVec> transferD_1P(4, dg::evaluate( dg::zero, grid));
    std::array<std::vector<double>, 4> transfer_prb;
    int dim_prb_ids[2] = {EtimeID, probeID};
    size_t count_prb[] = {1, probes.size()};
    size_t start_prb[] = {0, 0};
    int dataIDs_prb[4];
    std::string names_prb[4] = {"electrons_probes", "ions_probes", "potential_probes", "vorticity_probes"};

    if (p.save_pb){
        //Fields center point

        //////////////first output ////////////
        for (auto probe: probes){
            transfer_prb[0].push_back(y0[0][probe]);
            transfer_prb[1].push_back(y0[0][probe]);
            transfer_prb[2].push_back(exp.potential()[probe]);
            transfer_prb[3].push_back(y0[1][probe]);
        }

        for( unsigned i=0; i<4; i++){
            err = nc_def_var( ncid, names_prb[i].data(),  NC_DOUBLE, 2,  dim_prb_ids, &dataIDs_prb[i]);
            err = nc_put_vara_double( ncid, dataIDs_prb[i], start_prb, count_prb, transfer_prb[i].data());
        }

        err = nc_close(ncid);
    }

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

            Estart[0]    += 1;
            // start[0] +=1;
            if (p.save_pb){start_prb[0] += 1;}
            {
                err = nc_open(argv[2], NC_WRITE, &ncid);
                err = nc_put_vara_double( ncid, EtimevarID, Estart, Ecount, &time);
                if (p.save_pb)
                {
                    for (unsigned i = 0; i < probes.size(); i++){
                        transfer_prb[0][i] = y0[0][probes[i]];
                        transfer_prb[1][i] = y0[0][probes[i]];
                        transfer_prb[2][i] = exp.potential()[probes[i]];
                        transfer_prb[3][i] = y0[1][probes[i]];}
                    for (unsigned i = 0; i < 4; i++){
                        err = nc_put_vara_double( ncid, dataIDs_prb[i], start_prb, count_prb, transfer_prb[i].data());
                    }
                }

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
