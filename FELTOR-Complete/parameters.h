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
    std::string Time_Step;
	unsigned n, Nx, Ny;
	double dt;
	unsigned n_out, Nx_out, Ny_out;
	unsigned itstp;
	unsigned maxout;
	bool save_pb;
	std::vector<std::array<unsigned, 2>> probes;
	unsigned stages;

	double eps_pol, eps_time;
	double kappa, g, alpha, sgm, nu;
	double tau, nb, lambda;

	double tanh_width, profamp;
	double amp,  sigma,  posX,  posY;
	double amp2, sigma2, posX2, posY2;

	double lx, ly;
	double border_width, x_a, x_b, x_c;
	dg::bc bc_x_n, bc_y_n;
	dg::bc bc_x_omega, bc_y_omega;
	dg::bc bc_x_phi, bc_y_phi;


	Parameters( const Json::Value& js) {
	model     = js["model"].asString();
	modified  = js["modified"].asBool();
	save_pb   = js["save_probes"].asBool();
    Time_Step = js["Time_Step"].asString();
	n       = js["n"].asUInt();
	Nx      = js["Nx"].asUInt();
	Ny      = js["Ny"].asUInt();
	dt      = js["dt"].asDouble();
	n_out   = js["n_out"].asUInt();
	Nx_out  = js["Nx_out"].asUInt();
	Ny_out  = js["Ny_out"].asUInt();
	itstp   = js["itstp"].asUInt();
	maxout  = js["maxout"].asUInt();

	if (save_pb){
	try{
		for (auto probe: js["probes"])
			probes.push_back({probe[0].asUInt(), probe[1].asUInt()});
	}
	catch(...){
		probes.push_back({js["probes"][0].asUInt(), js["probes"][1].asUInt()});
	}}
	eps_pol  = js["eps_pol"].asDouble();
	eps_time = js["eps_time"].asDouble();
	stages   = js.get("stages",3).asUInt();
	kappa    = js["curvature"].asDouble();
	g        = js["dens_prof"].asDouble();
	alpha    = js["adiabatic"].asDouble();
	sgm      = js["sheath_diss"].asDouble();
	lambda   = js["sheath_pot"].asDouble();
	tau      = js["tau"].asDouble();
	nb       = js["nb"].asDouble();
	nu       = js["nu_perp"].asDouble();     // Dissipation

	tanh_width = js["tanh_width"].asDouble();
	profamp    = js["ProfAmp"].asDouble();
	amp        = js["amplitude"].asDouble();
	sigma      = js["sigma"].asDouble();
	posX       = js["posX"].asDouble();
	posY       = js["posY"].asDouble();

	lx           = js["lx"].asDouble();
	ly           = js["ly"].asDouble();
	border_width = js["border_width"].asDouble();
	x_a          = js["x_a"].asDouble() * lx;
	x_b          = js["x_b"].asDouble() * lx;
	x_c          = js["x_c"].asDouble() * lx;

	bc_x_n      = dg::str2bc(js["bc_x_n"].asString());
	bc_y_n      = dg::str2bc(js["bc_y_n"].asString());
	bc_x_omega  = dg::str2bc(js["bc_x_omega"].asString());
	bc_y_omega  = dg::str2bc(js["bc_y_omega"].asString());
	bc_x_phi    = dg::str2bc(js["bc_x_phi"].asString());
	bc_y_phi    = dg::str2bc(js["bc_y_phi"].asString());
	}

    void display( std::ostream& os = std::cout ) const
    {   os << "The model we are using is " <<model<<"\n"
			<< "Use the modified HW model: " << modified<<"\n"
            << "    With time stepper:    "<< Time_Step <<"\n";

		os << "Position of the boundaries between models \n"
			<< "         x_a:          = " << x_a << "\n"
			<< "         x_b:          = " << x_b << "\n"
			<< "         x_c:          = " << x_c << "\n"
			<< "width of the boundary: = " << border_width << "\n";

        os << "Physical parameters are: \n"
            <<"              Viscosity:              = "<<nu<<"\n"
            <<"              Curvature:              = "<<kappa<<"\n"
			<<"Density profile (not in use anymore): = "<<g<<"\n"
			<<"           adiabatic paramete:        = "<<alpha<<"\n"
			<<"           Sheath dissipation:        = "<<sgm<<"\n"
			<<"            Sheath potential:         = "<<lambda<<"\n"
			<<"             Source's tau             = "<<tau<<"\n"
			<<"                 nb                   = "<<nb<<"\n";

		if (save_pb){
		for(auto probe: probes){
			os << "Position measured at " << probe[0]
			   << "and " << probe[1] << "\n";}}

		os  << "Width of the initial profile: " << tanh_width << "\n"
			<< "        With amplitude:       " << profamp << "\n";

        os  <<"Blob parameters are: \n"
            << "     width:     = "<<sigma<<"\n"
            << "   amplitude:   = "<<amp<<"\n"
            << "     posX:      = "<<posX<<"\n"
            << "     posY:      = "<<posY<<"\n";

        os << "Boundary parameters are: \n"
            <<"    lx = "<<lx<<"\n"
            <<"    ly = "<<ly<<"\n";

			 //Curious! dg:: is not needed due to ADL!
        os << "Boundary conditions in x are: \n"
			<<"   n:  "<<bc2str(bc_x_n)<<"\n"
			<<"omega: "<<bc2str(bc_x_omega)<<"\n"
			<<"  phi: "<<bc2str(bc_x_phi)<<"\n";
        os << "Boundary conditions in y are: \n"
			<<"   n:  "<<bc2str(bc_y_n)<<"\n"
			<<"omega: "<<bc2str(bc_y_omega)<<"\n"
			<<"  phi: "<<bc2str(bc_y_phi)<<"\n";

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
