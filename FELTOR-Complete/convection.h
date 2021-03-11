#pragma once
#include <exception>

#include "dg/algorithm.h"
#include "parameters.h"

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
    const double m_nu, m_kappa;
    const double m_g, m_fau;
    container m_alpha, m_sigma;
    container m_lambda, m_lmbd_phi, m_nb;
    const container m_x, m_vol2d;

    container m_phi, m_temp, m_phi_perturbation, m_n_perturbation;
    container m_dn_dy;
    container m_exp_phi;
    container m_vx, m_vy;
    std::array<container,2> m_lapy;

    //matrices and solvers
    Matrix m_dy_n, m_dy_phi;
    Matrix         m_dx_phi;
    dg::Advection <Geometry, Matrix, container> m_advection_n, m_advection_omega;
    dg::Elliptic<Geometry, Matrix, container> m_laplaceM; //contains negative laplacian

    dg::MultigridCG2d<Geometry, Matrix, container> m_multigrid;
    std::vector<dg::Elliptic<Geometry, Matrix, container> > m_multi_pol;
    dg::Extrapolation<container> m_old_phi;

    dg::Average<container> m_average;
    const container m_S_domain;
    container m_Source;
    std::array<double,4> m_invariant, m_invariant_diss;
};

template< class Geometry, class M, class container>
ExplicitPart< Geometry, M, container>::ExplicitPart( const Geometry& grid, const Parameters& p ):
    m_model(p.modified), m_modified(p.modified),
    m_eps_pol(p.eps_pol),  m_nu(p.nu), m_kappa(p.kappa), m_g(p.g),
    m_fau(1. / p.tau),
    m_alpha(dg::evaluate(dg::one, grid)),
    m_sigma(m_alpha),
    m_lambda(m_alpha), m_lmbd_phi(m_alpha), m_nb(m_alpha),
    m_x( dg::evaluate( dg::cooX2d, grid)), m_vol2d( dg::create::volume(grid)),
    m_phi( evaluate( dg::zero, grid)), m_temp(m_phi), m_phi_perturbation(m_phi),
    m_n_perturbation(m_phi), m_dn_dy(m_phi), m_exp_phi(m_phi),
    m_vx(m_phi), m_vy(m_phi),
    m_lapy({ m_phi, m_phi}),
    m_dy_n( dg::create::dy(grid,   p.bc_y_n)),
    m_dy_phi( dg::create::dy(grid, p.bc_y_phi, dg::centered)),
    m_dx_phi( dg::create::dy(grid, p.bc_x_phi, dg::centered)),
    m_advection_n(     grid, p.bc_x_n,    p.bc_y_n),
    m_advection_omega( grid, p.bc_x_omega, p.bc_y_omega),
    m_laplaceM( grid, dg::normed, dg::centered),
    m_multigrid( grid, p.stages),
    m_old_phi( 2, m_phi),
    m_average(grid, dg::coo2d::y),
    m_S_domain(dg::evaluate(dg::TanhProfX(p.x_a, 0.001, -1., 0., 1.), grid)),
    m_Source(m_S_domain)
    {
    // m_g += -p.kappa;
    //construct multigrid
    m_multi_pol.resize(p.stages);
    for( unsigned u=0; u<p.stages; u++)
        m_multi_pol[u].construct( m_multigrid.grid(u), dg::not_normed, dg::centered);

    if (p.model == "IC_HW"){
      m_alpha = dg::evaluate(dg::TanhProfX(p.x_b, p.tanh_width,  1., 0., 1.), grid);
      m_sigma = dg::evaluate(dg::TanhProfX(p.x_b, p.tanh_width, -1., 0., 1.), grid);
    }
    if (p.model == "HW_IC"){
      m_alpha = dg::evaluate(dg::TanhProfX(p.x_b, p.tanh_width, -1., 0., 1.), grid);
      m_sigma = dg::evaluate(dg::TanhProfX(p.x_b, p.tanh_width,  1., 0., 1.), grid);
    }
    dg::blas1::scal(m_alpha,  p.alpha);
    dg::blas1::scal(m_sigma,  p.sgm);
    dg::blas1::scal(m_lambda, p.lambda);
    dg::blas1::scal(m_nb,     p.nb * m_fau);

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
    if(  number[0] == m_multigrid.max_iter())
		//// This gives an accuracy not reached fail
        throw dg::Fail( m_eps_pol);

	/// insert phi for the next extrapolation
    m_old_phi.update( m_phi);

    // v_x  = -dy phi (phi is defined negative)
    dg::blas2::symv( 1., m_dy_phi, m_phi, 0., m_vx);
    // v_y = dx phi (phi is defined negative)
    dg::blas2::symv(-1., m_dx_phi, m_phi, 0., m_vy);

    // Create the Exponential
    dg::blas1::axpby(1., m_lambda, 1., m_phi, m_lmbd_phi);
    dg::blas1::transform( m_lmbd_phi, m_exp_phi, dg::EXP<double>());
    // m_exp_phi =  dg::evaluate(dg::ExpProfX(1., 0., -1.), m_lmbd_phi);
    dg::blas1::pointwiseDivide(m_exp_phi, m_sigma, m_exp_phi);


	/////////////////////////update energetics/////////////////////
    for( unsigned i=0; i<2; i++)
	/// Calculate the Laplacian of y and place it in m_lapy
        ///              M         * x   = y
	dg::blas2::symv( m_laplaceM, y[i], m_lapy[i]);

	////  M = M * a
    dg::blas1::scal( m_lapy, -1.);

	//mass inveriant

	/// Total mass: Integrate n in the V    Vol       n
    m_invariant[0]      =  dg::blas1::dot( m_vol2d, y[0] );
	/// Diffusion of the total mass
    m_invariant_diss[0] = m_nu*dg::blas1::dot( m_vol2d, m_lapy[0]);

	//energy terms
	/// Total entropy
    m_invariant[1] = 0.5*dg::blas2::dot( y[0], m_vol2d, y[0]);

	//// Calculate Temperature, I believe
    m_laplaceM.variation( m_phi, m_temp);
	/// Total Kinetic energy, associated with the thermal energy
	/// the contribution is 1/2 for each degreed of freedom 1 ?? (3 ??), in this case
    m_invariant[2] = 0.5*dg::blas1::dot( m_vol2d, m_temp);

	/// Total potential energy, -Kappa * X coor  *  Volum  * n
    m_invariant[3] = -m_kappa * dg::blas2::dot( m_x, m_vol2d, y[0]);

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
        m_advection_n.upwind(     -1., m_vx, m_vy, y[0], 0., yp[0]);
        m_advection_omega.upwind( -1., m_vx, m_vy, y[1], 0., yp[1]);
        for (unsigned i=0; i<2;i++) {
            dg::blas1::pointwiseDot ( -1., m_alpha,   m_n_perturbation, 1., yp[i]);
            dg::blas1::pointwiseDot ( -1., m_alpha, m_phi_perturbation, 1., yp[i]);
    }}
	else {
            ///// [m_phi, y] = yp We compute the Poisson brackets
        m_advection_n.upwind(     -1., m_vx, m_vy, y[0], 0., yp[0]);
        m_advection_omega.upwind( -1., m_vx, m_vy, y[1], 0., yp[1]);
	}
    /// Kappa * d / dy, the term of the y derivative,
	/// phi is negative, that's why both have same sign
  /// IMPORTANT: kappa is not a salar anymore
	///              -m_kappa * M   *   x + a * y
    dg::blas2::gemv(m_dy_n, y[0], m_dn_dy);   // n derivative
    dg::blas1::axpby(-m_g, m_dn_dy, 1., yp[0]);
                               // -g *  n  *  d_y phi (V_x)
    dg::blas1::pointwiseDot(    -m_g,  y[0], m_vx, 1., yp[0]);
                              // g   d_y n  /  n
    dg::blas1::pointwiseDivide( -m_g, m_dn_dy, y[0], 1., yp[1]);
    // Exponentials
    dg::blas1::axpby(-1, m_exp_phi, 1., yp[0]);
    dg::blas1::axpbypgz(1, m_sigma, -1, m_exp_phi, 1., yp[1]);
    // Source
    dg::blas1::axpbypgz(m_fau, y[0], -1., m_nb, 0., m_Source);
    dg::blas1::pointwiseDot(-1., m_S_domain, m_Source, 1., yp[0]);

    return;
}

}//namespace convection
