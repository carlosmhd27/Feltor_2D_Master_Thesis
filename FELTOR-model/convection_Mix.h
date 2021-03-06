#pragma once
#include <exception>
// #include <string.h>
// #include <stdio.h>

#include "dg/algorithm.h"
#include "parameters_Mix.h"

namespace convection
{

template<class Geometry, class Matrix, class container>
struct ImplicitPart
{
    ImplicitPart( const Geometry& g, double nu):
        m_nu(nu),
        m_LaplacianM_perp( g, dg::normed, dg::centered),
        m_temp( evaluate( dg::zero, g))
    { }
    void operator()(double t, const std::array<container,2>& y, std::array<container,2>& yp)
    {
        // y[0] := n
        // y[1] := omega
        for( unsigned i=0; i<2; i++)
        {
            dg::blas2::symv( m_LaplacianM_perp, y[i], m_temp);
            dg::blas2::symv( m_LaplacianM_perp, m_temp, yp[i]);
        }
        dg::blas1::scal( yp, -m_nu);
    }

    const container& weights() const {
        return m_LaplacianM_perp.weights();
    }
    const container& inv_weights() const {
        return m_LaplacianM_perp.inv_weights();
    }
    const container& precond(){
        return m_LaplacianM_perp.precond();
    }

  private:
    double m_nu;
    dg::Elliptic<Geometry, Matrix, container> m_LaplacianM_perp;
    container m_temp;
};

template< class Geometry,  class Matrix, class container >
struct ExplicitPart
{
    ExplicitPart( const Geometry& g, const Parameters& p );

    const container& potential( ) const {
        return m_phi;
    }

    void operator()( double t, const std::array<container,2>& y, std::array<container,2>& yp);

    std::array<double,4> invariants( ) const {
        //0 : mass
        //1 : entropy
        //2 : kinetic
        //3 : potential
        return m_invariant;
    }
    std::array<double, 4> invariants_diffusion( ) const{
        //0 : mass
        //1 : entropy
        //2 : kinetic
        //3 : potential
        return m_invariant_diss;
    }

  private:
    bool  m_model;
    const bool m_modified;
    const double m_eps_pol;
    const double m_nu, m_kappa_var;
    double m_g;
    container m_alpha, m_kappa;
    const container m_x, m_vol2d;

    container m_phi, m_temp, m_phi_perturbation, m_n_perturbation;
    container m_dy_y0, m_dy_phi;
    std::array<container,2> m_lapy;

    //matrices and solvers
    Matrix m_dy;
    dg::ArakawaX< Geometry, Matrix, container> m_arakawa;
    dg::Elliptic<Geometry, Matrix, container> m_laplaceM; //contains negative laplacian

    dg::MultigridCG2d<Geometry, Matrix, container> m_multigrid;
    std::vector<dg::Elliptic<Geometry, Matrix, container> > m_multi_pol;
    dg::Extrapolation<container> m_old_phi;

    dg::Average<container> m_average;

    std::array<double,4> m_invariant, m_invariant_diss;
};

template< class Geometry, class M, class container>
ExplicitPart< Geometry, M, container>::ExplicitPart( const Geometry& grid, const Parameters& p ):
    m_model(p.modified), m_modified(p.modified),
    m_eps_pol(p.eps_pol),  m_nu(p.nu), m_kappa_var(p.kappa), m_g(p.g),
    m_alpha(dg::evaluate(dg::one, grid)),
    m_kappa(dg::evaluate(dg::one, grid)),
    m_x( dg::evaluate( dg::cooX2d, grid)), m_vol2d( dg::create::volume(grid)),
    m_phi( evaluate( dg::zero, grid)), m_temp(m_phi), m_phi_perturbation(m_phi),
    m_n_perturbation(m_phi), m_dy_y0(m_phi), m_dy_phi(m_phi),
    m_lapy({ m_phi, m_phi}),
    m_dy( dg::create::dy(grid)),
    m_arakawa( grid),
    m_laplaceM( grid, dg::normed, dg::centered),
    m_multigrid( grid, p.stages),
    m_old_phi( 2, m_phi),
    m_average(grid, dg::coo2d::y)
{
    m_g += -p.kappa;
    //construct multigrid
    m_multi_pol.resize(p.stages);
    for( unsigned u=0; u<p.stages; u++)
        m_multi_pol[u].construct( m_multigrid.grid(u), dg::not_normed, dg::centered);

    // IMPORTANT: Verify this and find a way to give width and sign
    // Also, dont use m_model, new variable as Mix = 'tanh' or 'step' or 'No'
    if (p.model == "IC_HW"){
      m_alpha = dg::evaluate(dg::TanhProfX(p.x_b, p.tanh_width,  1., 0., 1.), grid);
      m_kappa = dg::evaluate(dg::TanhProfX(p.x_b, p.tanh_width, -1., 0., 1.), grid);
    }
    if (p.model == "HW_IC"){
      m_alpha = dg::evaluate(dg::TanhProfX(p.x_b, p.tanh_width, -1., 0., 1.), grid);
      m_kappa = dg::evaluate(dg::TanhProfX(p.x_b, p.tanh_width,  1., 0., 1.), grid);
    }
    dg::blas1::scal(m_alpha, p.alpha);
    dg::blas1::scal(m_kappa, p.kappa);

    if (p.model.find("HW") != std::string::npos){m_model = true;}
    else {m_model = false;}
}

template< class G, class M, class container>
void ExplicitPart<G, M, container>::operator()( double t, const std::array<container,2>& y, std::array<container,2>& yp)
{
    //y[0] == n
    //y[1] == omega

	//yp[0] == [phi, n]
	//yp[1] == [phi, omega]

    /////////////////First, invert polarisation equation///////////
    //Note that we get the negative potential!!
	/// m_old_phi is an extrapolation class
    m_old_phi.extrapolate( m_phi);

	/// m_multigrid is a Conjugate Gradient in 2D for the multigrid
	/// m_multi_pol seems to be the Laplacian for the multigrid
	/// This is used to get the number of iterations needed to obtain
	/// the desired accuracy, the operation used for this is   operation (phi) = omega * weights
	/// where weights must be in the grid
    std::vector<unsigned> number = m_multigrid.direct_solve( m_multi_pol, m_phi, y[1], {m_eps_pol, 10*m_eps_pol, 10*m_eps_pol});

	/// insert phi for the next extrapolation
    m_old_phi.update( m_phi);
    if(  number[0] == m_multigrid.max_iter())
		//// This gives an accuracy not reached fail
        throw dg::Fail( m_eps_pol);

	/////////////////////////update energetics/////////////////////
  for( unsigned i=0; i<2; i++)
  	/// Calculate the Laplacian of y and place it in m_lapy
          ///              M         * x   = y
  	dg::blas2::symv( m_laplaceM, y[i], m_lapy[i]);

		////  M = M * a
    dg::blas1::scal( m_lapy, -1.);

	//mass inveriant

	/// Total mass: Integrate n in the V    Vol       n
    m_invariant[0]      =        dg::blas1::dot( m_vol2d, y[0] );
	/// Diffusion of the total mass
    m_invariant_diss[0] = m_nu * dg::blas1::dot( m_vol2d, m_lapy[0]);

	//energy terms
	/// Total entropy
    m_invariant[1] = 0.5*dg::blas2::dot( y[0], m_vol2d, y[0]);

	//// Calculate Temperature, I believe
    m_laplaceM.variation( m_phi, m_temp);
	/// Total Kinetic energy, associated with the thermal energy
	/// the contribution is 1/2 for each degreed of freedom 1 ?? (3 ??), in this case
    m_invariant[2] = 0.5*dg::blas1::dot( m_vol2d, m_temp);

    /// IMPORTANT: kappa is not a salar anymore
	  /// Total potential energy, -Kappa * X coor  *  Volum  * n
    m_invariant[3] = - m_kappa_var * dg::blas2::dot( m_x, m_vol2d, y[0]);

    //energy dissipation terms, the same over the Laplacian of y
    m_invariant_diss[0] =  m_nu*dg::blas1::dot( m_vol2d, m_lapy[0]);
    m_invariant_diss[1] =  m_nu*dg::blas2::dot( y[0],   m_vol2d, m_lapy[0]);
    m_invariant_diss[2] = -m_nu*dg::blas2::dot( m_phi, m_vol2d, m_lapy[1]);
    m_invariant_diss[3] = -m_nu*dg::blas2::dot( m_x,   m_vol2d, m_lapy[0]);
    ///////////////////////Equations////////////////////////////////
	if (m_model) {
		///Average///
	    if(m_modified){
            m_average( m_phi, m_phi_perturbation);
            m_average(  y[0], m_n_perturbation);
            dg::blas1::axpby( 1., m_phi, -1., m_phi_perturbation);
            dg::blas1::axpby( 1.,  y[0], -1., m_n_perturbation);
        }
        else
        {
            dg::blas1::copy( m_phi, m_phi_perturbation);
            dg::blas1::copy( y[0], m_n_perturbation);
        }
        ///Perturbation terms///
        for (unsigned i=0; i<2;i++) {
            m_arakawa( m_phi, y[i], yp[i]);
            dg::blas1::pointwiseDot ( -1., m_alpha,   m_n_perturbation, 1., yp[i]);
            dg::blas1::pointwiseDot ( -1., m_alpha, m_phi_perturbation, 1., yp[i]);
    }}
	else {
    	for( unsigned i=0; i<2; i++)
			///// [m_phi, y] = yp
    	    m_arakawa( m_phi,y[i], yp[i]); //m_phi is negative!
	}
	/// Kappa * d / dy, the term of the y derivative,
	/// phi is negative, that's why both have same sign
  /// IMPORTANT: kappa is not a salar anymore
	///              -m_kappa * M   *   x + a * y
    dg::blas2::gemv(m_dy, y[0], m_dy_phi);
    dg::blas2::gemv(m_dy, y[0], m_dy_y0);
    dg::blas2::axpby(m_g, m_dy_phi, 1., yp[0]);
    dg::blas2::pointwiseDot( -1., m_kappa, m_dy_phi, 1., yp[0]);
    dg::blas1::pointwiseDot( -1., m_kappa, m_dy_y0, 1., yp[1]);
    return;
}

}//namespace convection
