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
        m_LaplacianM_perp( g, dg::normed, dg::centered)
    { }
    void operator()(double t, const std::array<container,2>& y, std::array<container,2>& yp)
    {
        // y[0] := n
        // y[1] := omega
        for( unsigned i=0; i<2; i++)
            dg::blas2::symv( m_LaplacianM_perp, y[i], yp[i]);
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
    const double m_eps_pol;
    const double m_kappa, m_nu;
    const container m_x, m_vol2d;

    container m_phi, m_temp;
    std::array<container,2> m_lapy;

    //matrices and solvers
    Matrix m_dy;
    dg::ArakawaX< Geometry, Matrix, container> m_arakawa;
    dg::Elliptic<Geometry, Matrix, container> m_laplaceM; //contains negative laplacian

    dg::MultigridCG2d<Geometry, Matrix, container> m_multigrid;
    std::vector<dg::Elliptic<Geometry, Matrix, container> > m_multi_pol;
    dg::Extrapolation<container> m_old_phi;

    std::array<double,4> m_invariant, m_invariant_diss;

};

template< class Geometry, class M, class container>
ExplicitPart< Geometry, M, container>::ExplicitPart( const Geometry& grid, const Parameters& p ):
    m_eps_pol(p.eps_pol), m_kappa(p.kappa), m_nu(p.nu),
    m_x( dg::evaluate( dg::cooX2d, grid)), m_vol2d( dg::create::volume(grid)),
    m_phi( evaluate( dg::zero, grid)), m_temp(m_phi),
    m_lapy({ m_phi, m_phi}),
    m_dy( dg::create::dy(grid)),
    m_arakawa( grid),
    m_laplaceM( grid, dg::normed, dg::centered),
    m_multigrid( grid, p.stages),
    m_old_phi( 2, m_phi)
{
    //construct multigrid
    m_multi_pol.resize(p.stages);
    for( unsigned u=0; u<p.stages; u++)
        m_multi_pol[u].construct( m_multigrid.grid(u), dg::not_normed, dg::centered);
}

template< class G, class M, class container>
void ExplicitPart<G, M, container>::operator()( double t, const std::array<container,2>& y, std::array<container,2>& yp)
{
    //y[0] == n
    //y[1] == omega

    /////////////////First, invert polarisation equation///////////
    //Note that we get the negative potential!!
    m_old_phi.extrapolate( m_phi);
    std::vector<unsigned> number = m_multigrid.direct_solve( m_multi_pol, m_phi, y[1], m_eps_pol);
    m_old_phi.update( m_phi);
    if(  number[0] == m_multigrid.max_iter())
        throw dg::Fail( m_eps_pol);
    /////////////////////////update energetics/////////////////////
    for( unsigned i=0; i<2; i++)
        dg::blas2::symv( m_laplaceM, y[i], m_lapy[i]);
    dg::blas1::scal( m_lapy, -1.);
    //mass inveriant
    m_invariant[0]      =  dg::blas1::dot( m_vol2d, y[0] );
    m_invariant_diss[0] = m_nu*dg::blas1::dot( m_vol2d, m_lapy[0]);
    //energy terms
    m_invariant[1] = 0.5*dg::blas2::dot( y[0], m_vol2d, y[0]);
    m_arakawa.variation( m_phi, m_temp);
    m_invariant[2] = 0.5*dg::blas1::dot( m_vol2d, m_temp);
    m_invariant[3] = -m_kappa*dg::blas2::dot( m_x, m_vol2d, y[0]);
    //energy dissipation terms
    m_invariant_diss[0] =  m_nu*dg::blas1::dot( m_vol2d, m_lapy[0]);
    m_invariant_diss[1] =  m_nu*dg::blas2::dot( y[0],   m_vol2d, m_lapy[0]);
    m_invariant_diss[2] = -m_nu*dg::blas2::dot( m_phi, m_vol2d, m_lapy[1]);
    m_invariant_diss[3] = -m_nu*dg::blas2::dot( m_x,   m_vol2d, m_lapy[0]);
    ///////////////////////Equations////////////////////////////////
    for( unsigned i=0; i<2; i++)
        m_arakawa( m_phi,y[i], yp[i]); //m_phi is negative!

    dg::blas2::gemv( -m_kappa, m_dy, m_phi, 1., yp[0]);
    dg::blas2::gemv( -m_kappa, m_dy,  y[0], 1., yp[1]);
    return;
}

}//namespace convection
