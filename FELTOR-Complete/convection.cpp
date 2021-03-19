#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>

#include "draw/host_window.h"
//#include "draw/device_window.cuh"

#include "dg/algorithm.h"
#include "convection.h"
#include "parameters.h"

/*
   - reads parameters from input.json or any other given file,
   - integrates the convection - functor and
   - directly visualizes results on the screen
*/


int main( int argc, char* argv[])
{
    ////Parameter initialisation ////////////////////////////////////////////
    std::stringstream title;
    Json::Value js;
    if( argc == 1)
    {
        std::ifstream is("input.json");
        is >> js;
    }
    else if( argc == 2)
    {
        std::ifstream is(argv[1]);
        is >> js;
    }
    else
    {
        std::cerr << "ERROR: Too many arguments!\nUsage: "<< argv[0]<<" [filename]\n";
        return -1;
    }
    const convection::Parameters p( js);
    p.display( std::cout);
    /////////glfw initialisation ////////////////////////////////////////////
    GLFWwindow* w = draw::glfwInitAndCreateWindow( js.get("width",500.).asDouble(), js.get("height",1000).asDouble(), "");
    draw::RenderHostData render(js.get("rows",2).asDouble(), js.get("cols",1).asDouble());
    /////////////////////////////////////////////////////////////////////////
    dg::Grid2d grid( 0, p.lx, 0, p.ly, p.n, p.Nx, p.Ny, p.bc_x_n, p.bc_y_n);
    //create RHS
    convection::ExplicitPart<dg::CartesianGrid2d, dg::DMatrix, dg::DVec> exp( grid, p);
    convection::ImplicitPart<dg::CartesianGrid2d, dg::DMatrix, dg::DVec> imp( grid, p);
    //////////////////create initial vector///////////////////////////////////////
    dg::Gaussian g( p.posX*p.lx, p.posY*p.ly, p.sigma, p.sigma, p.amp); //gaussian width is in absolute values
    std::array<dg::DVec,2> y0{
        dg::evaluate(g, grid),
        dg::evaluate(dg::zero,grid) // omega == 0
    };
    {dg::DVec ones(dg::evaluate(dg::one, grid));
    dg::blas1::axpby( p.nb,  ones, 1., y0[0]);}
    //////////////////////////////////////////////////////////////////////
    dg::Karniadakis< std::array<dg::DVec,2> > stepper( y0, y0[0].size(), p.eps_time);

    dg::DVec dvisual( grid.size(), 0.);
    dg::HVec hvisual( grid.size(), 0.), visual(hvisual);
    dg::IHMatrix equi = dg::create::backscatter( grid);
    draw::ColorMapRedBlueExt colors( 1.);
    //create timer
    dg::Timer t;
    double time = 0;
    stepper.init( exp, imp, time, y0, p.dt);
    const double mass0 = exp.mass(), mass_blob0 = mass0 - grid.lx()*grid.ly();
    // double E0 = exp.invariants()[1] + exp.invariants()[2], energy0 = E0, E1 = 0, diff = 0;
    std::cout << "Begin computation \n";
    std::cout << std::scientific << std::setprecision( 2);
    unsigned step = 0;
    while ( !glfwWindowShouldClose( w ))
    {
        //transform n to an equidistant grid
        dvisual = y0[0];
        dg::assign( dvisual, hvisual);
        dg::blas2::gemv( equi, hvisual, visual);
        //compute the color scale
        colors.scale() =  (float)thrust::reduce( visual.begin(), visual.end(), 0., dg::AbsMax<double>() );
        title << std::setprecision(2) << std::scientific;
        title <<"n / "<<colors.scale()<<"\t";
        //draw and swap buffers
        render.renderQuad( visual, grid.n()*grid.Nx(), grid.n()*grid.Ny(), colors);

        //transform omega to an equidistant grid
        dg::assign( y0[1], hvisual);
        dg::blas2::gemv( equi, hvisual, visual);
        //compute the color scale
        colors.scale() =  (float)thrust::reduce( visual.begin(), visual.end(), 0., dg::AbsMax<double>() );
        //draw and swap buffers
        title <<"omega / "<<colors.scale()<<"\t";
        title << std::fixed;
        title << " &&   time = "<<time;
        render.renderQuad( visual, grid.n()*grid.Nx(), grid.n()*grid.Ny(), colors);
        glfwSetWindowTitle(w,title.str().c_str());
        title.str("");
        glfwPollEvents();
        glfwSwapBuffers( w);

        //step
#ifdef DG_BENCHMARK
        t.tic();
#endif//DG_BENCHMARK
        for( unsigned i=0; i<p.itstp; i++)
        {
            step++;
            {
                std::cout << "(m_tot-m_0)/m_0: "<< (exp.mass()-mass0)/mass_blob0<<"\t";
            }
            try{ stepper.step( exp, imp, time, y0);}
            catch( dg::Fail& fail) {
                std::cerr << "CG failed to converge to "<<fail.epsilon()<<"\n";
                std::cerr << "Does Simulation respect CFL condition?\n";
                glfwSetWindowShouldClose( w, GL_TRUE);
                break;
            }
        }
#ifdef DG_BENCHMARK
        t.toc();
        std::cout << "\n\t Step "<<step;
        std::cout << "\n\t Average time for one step: "<<t.diff()/(double)p.itstp<<"s\n\n";
#endif//DG_BENCHMARK
    }
    glfwTerminate();
    ////////////////////////////////////////////////////////////////////

    return 0;

}
