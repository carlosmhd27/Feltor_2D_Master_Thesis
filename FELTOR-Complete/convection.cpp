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
    //////////////////create initial vector///////////////////////////////////////
    dg::Gaussian g ( p.posX *p.lx, p.posY * p.ly, p.sigma,  p.sigma,  p.amp ); //gaussian width is in absolute values
    // dg::Gaussian g2( p.posX2*p.lx, p.posY2 * p.ly, p.sigma2, p.sigma2, p.amp2); //gaussian width is in absolute values
    dg::TanhProfX tanh (p.x_b, p.tanh_width,  -1., 0.001 - p.nb, p.nb); // x_0, width, sign, B, A => B + A * 0.5(1 + sign tanh((x- x_0) / width))
    std::array<dg::DVec,2> y0{
        dg::evaluate(tanh, grid),
        dg::evaluate(dg::zero,grid) // omega == 0
    };
    { // dg::DVec ones(dg::evaluate(dg::one, grid));
    dg::DVec g2(dg::evaluate(g, grid));
    dg::blas1::axpby( 1., g2, 1., y0[0]);}
    //MW: our density variable is now y0 = n-nb
    //{dg::DVec ones(dg::evaluate(dg::one, grid));
    //dg::blas1::axpby( p.nb,  ones, 1., y0[0]);}
    //////////////////////////////////////////////////////////////////////
    //dg::ExplicitMultistep< std::array<dg::DVec,2> > stepper( "TVB-3-3", y0);
    dg::Adaptive<dg::ERKStep<std::array<dg::DVec, 2>>> stepper( "Bogacki-Shampine-4-2-3", y0);

    dg::DVec dvisual( grid.size(), 0.);
    dg::HVec hvisual( grid.size(), 0.), visual(hvisual);
    dg::IHMatrix equi = dg::create::backscatter( grid);
    draw::ColorMapRedBlueExt colors( 1.);
    //create timer
    dg::Timer t;
    double time = 0;
    //stepper.init( exp, time, y0, p.dt);
    const double mass0 = exp.mass(), mass_blob0 = mass0 - grid.lx()*grid.ly();
    // double E0 = exp.invariants()[1] + exp.invariants()[2], energy0 = E0, E1 = 0, diff = 0;
    std::cout << "Begin computation \n";
    std::cout << std::scientific << std::setprecision( 2);
    unsigned step = 0;
    double dt = p.dt;
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
            //try{ stepper.step( exp, time, y0);}
            try{ stepper.step( exp, time, y0, time, y0, dt, dg::pid_control, dg::l2norm, p.eps_time, p.eps_time);}
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
