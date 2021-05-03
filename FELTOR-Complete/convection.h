#pragma once
#include <exception>

#include "dg/algorithm.h"
#include "parameters.h"

namespace convection
{

template< class Geometry,  class Matrix, class container >
struct ExplicitPart
{
    ExplicitPart( const Geometry& g, const Parameters& p );

    const container& potential( ) const {
        return m_phi;
    }
    const container& vradial( ) const {
        return m_vx;
    }

    void operator()( double t, const std::array<container,2>& y, std::array<container,2>& yp);

    double mass( ) const {
        //0 : mass
        return m_mass;
    }
    private:

    Parameters m_p;
    bool  m_IC, m_HW;
    container m_alpha, m_sigma;
    container m_lambda, m_lmbd_phi;
    const container m_x, m_vol2d;

    container m_phi, m_temp, m_phi_perturbation, m_n_perturbation;
    container m_dn_dy, m_full_n;
    container m_exp_phi, m_S_domain, m_W_domain;
    container m_vx, m_vy;

    //matrices and solvers
    const Matrix m_dy_n, m_dy_phi;
    const Matrix         m_dx_phi;
    dg::Advection <Geometry, Matrix, container> m_advection_n, m_advection_omega;

    dg::MultigridCG2d<Geometry, Matrix, container> m_multigrid;
    std::vector<dg::Elliptic<Geometry, Matrix, container> > m_multi_pol;
    dg::Elliptic<Geometry, Matrix, container> m_laplaceM;
    dg::Extrapolation<container> m_old_phi;

    dg::Average<container> m_average;
    double m_mass;
};

template< class Geometry, class M, class container>
ExplicitPart< Geometry, M, container>::ExplicitPart( const Geometry& grid, const Parameters& p ):
    m_p(p),
    m_IC(false), m_HW(false),
    m_alpha(dg::evaluate(dg::zero, grid)),
    m_sigma(m_alpha),
    m_lambda(m_alpha), m_lmbd_phi(m_alpha),
    m_x( dg::evaluate( dg::cooX2d, grid)), m_vol2d( dg::create::volume(grid)),
    m_phi( evaluate( dg::zero, grid)), m_temp(m_phi), m_phi_perturbation(m_phi),
    m_n_perturbation(m_phi), m_dn_dy(m_phi), m_full_n(m_phi), m_exp_phi(m_phi),
    m_S_domain(dg::evaluate(dg::TanhProfX(p.x_a, p.border_width, -1., 0., 1.), grid)),
    m_W_domain(dg::evaluate(dg::TanhProfX(p.x_c, p.border_width, +1., 0., 1.), grid)),
    m_vx(m_phi), m_vy(m_phi),
    // m_lapy({ m_phi, m_phi}),
    m_dy_n(   dg::create::dy(grid, p.bc_y_n)),
    m_dy_phi( dg::create::dy(grid, p.bc_y_phi)), //dg::centered is by default
    m_dx_phi( dg::create::dx(grid, p.bc_x_phi)), //dg::centered is by default
    m_advection_n(     grid, p.bc_x_n,     p.bc_y_n),
    m_advection_omega( grid, p.bc_x_omega, p.bc_y_omega),
    m_multigrid( grid, p.stages),
    m_laplaceM( grid, dg::normed, dg::centered),
    m_old_phi( 2, m_phi),
    m_average(grid, dg::coo2d::y)
    {
    //construct multigrid
    m_multi_pol.resize(p.stages);
    for( unsigned u=0; u<p.stages; u++)
        m_multi_pol[u].construct( m_multigrid.grid(u), dg::not_normed, dg::centered);

    if (p.model == "IC_HW"){                               /// sign, B, A
      m_alpha = dg::evaluate(dg::TanhProfX(p.x_b, p.border_width,  1., 0., 1.), grid);
      m_sigma = dg::evaluate(dg::TanhProfX(p.x_b, p.border_width, -1., 0., 1.), grid);
    }
    if (p.model == "HW_IC"){
      m_alpha = dg::evaluate(dg::TanhProfX(p.x_b, p.border_width, -1., 0., 1.), grid);
      m_sigma = dg::evaluate(dg::TanhProfX(p.x_b, p.border_width,  1., 0., 1.), grid);
    }
    if (p.model == "IC"){                               /// sign, B, A
      m_alpha = dg::evaluate(dg::zero, grid);
      m_sigma = dg::evaluate(dg::one, grid);
    }
    if (p.model == "HW"){
      m_alpha = dg::evaluate(dg::one, grid);
      m_sigma = dg::evaluate(dg::zero, grid);
    }
    dg::blas1::scal(m_alpha,  p.alpha);
    dg::blas1::scal(m_sigma,  p.sgm);
    dg::blas1::scal(m_lambda, p.lambda);

    if (p.model.find("IC") != std::string::npos){m_IC = true;}
    else {m_IC = false;}

    if (p.model.find("HW") != std::string::npos){m_HW = true;}
    else {m_HW = false;}
    }

template< class G, class M, class container>
void ExplicitPart<G, M, container>::operator()( double t, const std::array<container,2>& y, std::array<container,2>& yp)
{
    //y[0] == n-nb
    //y[1] == omega

	//yp[0] == [dn/dt]
	//yp[1] == [domega/dt]

    /////////////////First, invert polarisation equation///////////
    //Note that we get the negative potential!!
	/// m_old_phi is an extrapolation class
    m_old_phi.extrapolate( m_phi);

	/// m_multigrid is a Conjugate Gradient in 2D for the multigrid
	/// m_multi_pol seems to be the Laplacian for the multigrid
	/// This is used to get the number of iterations needed to obtain
	/// the desired accuracy, the operation used for this is   operation (phi) = omega * weights
	/// where weights must be in the grid
    std::vector<unsigned> number = m_multigrid.direct_solve( m_multi_pol, m_phi, y[1], {m_p.eps_pol, m_p.eps_pol, m_p.eps_pol});
    if(  number[0] == m_multigrid.max_iter())
		//// This gives an accuracy not reached fail
        throw dg::Fail( m_p.eps_pol);

	/// insert phi for the next extrapolation
    m_old_phi.update( m_phi);
    dg::blas1::scal(m_phi, -1.);

    // v_x  = -dy phi (phi is defined negative)
    dg::blas2::symv( -1., m_dy_phi, m_phi, 0., m_vx);
    // v_y = dx phi (phi is defined negative)
    dg::blas2::symv(  1., m_dx_phi, m_phi, 0., m_vy);

	/// Total mass: Integrate n in the V    Vol       n
    m_mass =  dg::blas1::dot( m_vol2d, y[0] );

    ///////////////////////Equations////////////////////////////////
    if (m_HW) {
		///Average///
	    if(m_p.modified){
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
        //MW the total density advection includes curvature drift
        dg::blas1::transform(  m_vy, m_temp, dg::PLUS<double>(m_p.kappa));
        m_advection_n.upwind(     -1., m_vx, m_temp, y[0], 0., yp[0]);
        m_advection_omega.upwind( -1., m_vx, m_vy, y[1], 0., yp[1]);
        ///Perturbation terms///
        for (unsigned i=0; i<2;i++) {
            dg::blas1::pointwiseDot ( -1., m_alpha,   m_n_perturbation, 1., yp[i]);
            dg::blas1::pointwiseDot (  1., m_alpha, m_phi_perturbation, 1., yp[i]);
    }}
	else {
            ///// [m_phi, y] = yp We compute the Poisson brackets
            //MW the total density advection includes curvature drift
        dg::blas1::transform(  m_vy, m_temp, dg::PLUS<double>(m_p.kappa));
        m_advection_n.upwind(     -1., m_vx, m_temp, y[0], 0., yp[0]);
        m_advection_omega.upwind( -1., m_vx, m_vy, y[1], 0., yp[1]);
	}
    /// Kappa * d / dy, the term of the y derivative,
	/// phi is negative, that's why both have same sign
///////////////////////////Complete/////////////////////////////////////
    dg::blas1::transform( y[0], m_full_n, dg::PLUS<double>(m_p.nb));
	///              -m_p.kappa * M   *   x + a * y
    dg::blas2::symv(m_dy_n, y[0], m_dn_dy);   // n derivative
    dg::blas1::pointwiseDot( -m_p.kappa,  m_full_n, m_vx, 1., yp[0]); // -(g-kappa) n d_y phi (-V_x)
    dg::blas1::axpby( m_p.g,  m_vx, 1., yp[0]); // g d_y phi (-V_x)
    dg::blas1::pointwiseDivide( -m_p.kappa, m_dn_dy, m_full_n, 1., yp[1]);// -kappa   d_y n  /  n
    ///Exponentials
    ///Create the Exponential
    if(m_IC){
        dg::blas1::axpby(1., m_lambda, -1., m_phi, m_lmbd_phi);
        dg::blas1::transform( m_lmbd_phi, m_exp_phi, dg::EXP<double>());
        dg::blas1::pointwiseDot(m_exp_phi, m_sigma, m_exp_phi);
        dg::blas1::pointwiseDot(-1, m_full_n, m_exp_phi, 1., yp[0]);
        dg::blas1::axpbypgz(1, m_sigma, -1, m_exp_phi, 1., yp[1]);
    }
///////////////////////////Complete/////////////////////////////////////
/////////////////////////////OLD////////////////////////////////////////
    //dg::blas2::symv(m_dy_n, y[0], m_dn_dy);   // n derivative

    //dg::blas1::axpby( -m_p.kappa+m_p.g, m_vx,    1., yp[0]);

    //dg::blas1::axpby( -m_p.kappa, m_dn_dy, 1., yp[1]);

/////////////////////////////OLD////////////////////////////////////////
////////////////////MW Explicit ///////////////////////////
    // Source
    // better do it explicit since the tau values do not really become small
    // and m_nu should be small in all cases
    if (m_HW){
        //dg::blas1::axpbypgz(m_fau, y[0], -1., m_nb, 0., m_Source);
        //dg::blas1::pointwiseDot(-1., m_S_domain, m_Source, 1., yp[0]);
        dg::blas1::pointwiseDot( -1./m_p.tau, y[0], m_S_domain, 1., yp[0]);
        dg::blas1::pointwiseDot( -1./m_p.tau, y[1], m_W_domain, 1., yp[1]);
    }
    for( unsigned i=0; i<2; i++)
    {
        dg::blas2::symv( m_laplaceM, y[i], m_temp);
        dg::blas2::symv( -m_p.nu, m_laplaceM, m_temp, 1., yp[i]);
    }
    return;
}

}//namespace convection
