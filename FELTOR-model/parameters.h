#pragma once
#include <string>
#include "dg/algorithm.h"
// #include "json/json.h"
#include <jsoncpp/json/json.h>

namespace convection{

/**
 * @brief Provide a mapping between input file and named parameters
 */
struct Parameters
{
	std::string model;
	bool modified;
	unsigned n, Nx, Ny;
	double dt;
	unsigned n_out, Nx_out, Ny_out;
	unsigned itstp;
	unsigned maxout;
	unsigned save_pb;
	std::vector<std::array<unsigned, 2>> probes;
	unsigned stages;

	double eps_pol, eps_time;
	double kappa, g, alpha;
	std::array<double, 2> nu;

	double amp, sigma, posX, posY;

	double lx, ly;
	dg::bc bc_x, bc_y;

	Parameters( const Json::Value& js) {
	model    = js["model"].asString();
	modified = js["modified"].asBool();

	n       = js["n"].asUInt();
	Nx      = js["Nx"].asUInt();
	Ny      = js["Ny"].asUInt();
	dt      = js["dt"].asDouble();
	n_out   = js["n_out"].asUInt();
	Nx_out  = js["Nx_out"].asUInt();
	Ny_out  = js["Ny_out"].asUInt();
	itstp   = js["itstp"].asUInt();
	maxout  = js["maxout"].asUInt();
	save_pb = js["save_probes"].asUInt();

	if (save_pb){
	try{
		for (auto probe: js["probes"])
			probes.push_back({probe[0].asUInt(), probe[1].asUInt()});
	}
	catch(...){
		probes.push_back({js["probes"][0].asUInt(), js["probes"][1].asUInt()});
	}}
	eps_pol     = js["eps_pol"].asDouble();
	eps_time    = js["eps_time"].asDouble();
	stages      = js.get("stages",3).asUInt();
	kappa       = js["curvature"].asDouble();
	g           = js["dens_prof"].asDouble();
	alpha       = js["adiabatic"].asDouble();

	nu[0]       = js["nu_perp"].asDouble(); // Density dissipation
	if (js["mu_perp"].asDouble()){nu[1]= js["mu_perp"].asDouble();}
	else {nu[1]= js["nu_perp"].asDouble();} // Vorticity Dissipation

	amp         = js["amplitude"].asDouble();
	sigma       = js["sigma"].asDouble();
	posX        = js["posX"].asDouble();
	posY        = js["posY"].asDouble();
	lx          = js["lx"].asDouble();
	ly          = js["ly"].asDouble();
	bc_x        = dg::str2bc(js["bc_x"].asString());
	bc_y        = dg::str2bc(js["bc_y"].asString());
  }

    void display( std::ostream& os = std::cout ) const
    {   os << "The model we are using is " <<model<<"\n"
		    << "Use the modified HW model: " << modified<<"\n";
        os << "Physical parameters are: \n"
			<<" Density Viscosity:  = "<<nu[0]<<"\n"
			<<"Vorticity Viscosity: = "<<nu[1]<<"\n"
            <<"    Curvature:       = "<<kappa<<"\n"
			<<"  density profile:   = "<<g<<"\n"
			<<" adiabatic parameter = "<<alpha<<"\n";
        os  <<"Blob parameters are: \n"
            << "    width:        "<<sigma<<"\n"
            << "    amplitude:    "<<amp<<"\n"
            << "    posX:         "<<posX<<"\n"
            << "    posY:         "<<posY<<"\n";

		if (save_pb){
		for(auto probe: probes){
			os << "Position measured at " << probe[0]
			   << "and " << probe[1] << "\n";}}

        os << "Boundary parameters are: \n"
            <<"    lx = "<<lx<<"\n"
            <<"    ly = "<<ly<<"\n";
        os << "Boundary conditions in x are: \n"
            <<"    "<<bc2str(bc_x)<<"\n";  //Curious! dg:: is not needed due to ADL!
        os << "Boundary conditions in y are: \n"
            <<"    "<<bc2str(bc_y)<<"\n";
        os << "Algorithmic parameters are: \n"
            <<"    n  = "<<n<<"\n"
            <<"    Nx = "<<Nx<<"\n"
            <<"    Ny = "<<Ny<<"\n"
            <<"    dt = "<<dt<<"\n"
            << "   Accuracy for CG:          "<<eps_pol<<"\n"
            << "   Accuracy for Timestepper: "<<eps_time<<"\n"
            << "   # of multgrid stages:     "<<stages<<"\n";
        os << "Output parameters: \n"
            <<"    n_out  = "<<n_out<<"\n"
            <<"    Nx_out = "<<Nx_out<<"\n"
            <<"    Ny_out = "<<Ny_out<<"\n"
            <<"    Steps between outputs: "<<itstp<<"\n"
            <<"    Number of outputs:     "<<maxout<<std::endl; //the endl is for the implicit flush
    }
};
}//namespace convection
