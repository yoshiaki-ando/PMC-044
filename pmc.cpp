#include <iostream>
#include <complex>

#include <Mie_scattering.h>

#include "Fixed_parameter.h"
#include "pmc.h"

PMC::PMC(void){
  set_gauss_legendre();
}

PMC::PMC(
      double lambda, std::complex <double> m,
      AndoLab::Vector3d <double> rc, double sig_z, double sig_r,
      double r0, double sig,
      double N0):
          pR_center( rc ), pSigma_z( sig_z ), pSigma_r( sig_r ),
          pR0( r0 ), pSigma( sig ), pN0( N0 ),
          mie(lambda, r0, m)
        {
  set_gauss_legendre();
  calc_normalized_beta();
  pR_center_magnitude = pR_center.abs();
  dummy = false;
}

void PMC::set(double lambda, std::complex <double> m,
      AndoLab::Vector3d <double> rc, double sig_z, double sig_r,
      double r0, double sig,
      double N0){
  pR_center = rc;
  pSigma_z = sig_z;
  pSigma_r = sig_r;
  pR0 = r0;
  pSigma = sig;
  pN0 = N0;
  mie.set(lambda, r0, m);

  calc_normalized_beta();
  pR_center_magnitude = pR_center.abs();
  dummy = false;
}

void PMC::set_gauss_legendre(void){
  constexpr double zp[8] = {
      0.09501250983763744018531933542496,
      0.28160355077925891323046050146050,
      0.45801677765722738634241944298358,
      0.61787624440264374844667176404879,
      0.75540440835500303389510119484744,
      0.86563120238783174388046789771239,
      0.94457502307323257607798841553461,
      0.98940093499164993259615417345033 };
  constexpr double w[8] = {
      0.18945061045506849628539672320828,
      0.18260341504492358886676366796922,
      0.16915651939500253818931207903036,
      0.14959598881657673208150173054748,
      0.12462897125553387205247628219202,
      0.09515851168249278480992510760225,
      0.06225352393864789286284383699438,
      0.02715245941175409485178057245602 };

  for(int i = 0; i < 8; i++){
    gauss_z[i] = zp[i];
    gauss_w[i] = w[i];
  }
}

//void PMC::set_spatial_sd(double sig_r, double sig_z){
//  pSigma_r = sig_r;
//  pSigma_z = sig_z;
//}
//
/* ゼロで切断した粒径の正規分布 */
double truncated_normal
(
 const double r,
 const double r0,
 const double sigma
 ){

  const double K0 = 0.5 * ( 1.0 + std::erf(r0 / std::sqrt(2.0) / sigma) );
  const double x = (r - r0) / sigma;
  return std::exp(-0.5 * x * x) / std::sqrt(2.0 * M_PI) / sigma / K0;
}

/* 全波長にわたる散乱係数β (全数密度Nの係数は別) */
void PMC::calc_normalized_beta(void){

  double ans = 0.0;
  double ax = 0.5 * (radius_max - radius_min);
  double bx = 0.5 * (radius_max + radius_min);

  for(int i = 0; i < 8; i++){
    ans += gauss_w[i] * (
        truncated_normal( ax*gauss_z[i] + bx, pR0, pSigma )
        * mie.scs_t( ax*gauss_z[i] + bx)
        +
        truncated_normal( -ax*gauss_z[i] + bx, pR0, pSigma )
        * mie.scs_t( -ax*gauss_z[i] + bx)
        );
  }

  pNormalized_Beta = ans * ax;
}

/* 全波長にわたる散乱係数β(θ) (全数密度Nの係数は別) */
void PMC::calc_normalized_beta_th(double th){

  double ans = 0.0;
  double ax = 0.5 * (radius_max - radius_min);
  double bx = 0.5 * (radius_max + radius_min);

  for(int i = 0; i < 8; i++){
    ans += gauss_w[i] *
      (
       truncated_normal( ax*gauss_z[i] + bx, pR0, pSigma )
       * mie.scs( ax*gauss_z[i] + bx, th)
       +
       truncated_normal( -ax*gauss_z[i] + bx, pR0, pSigma )
       * mie.scs( -ax*gauss_z[i] + bx, th)
       );
  }

  pNormalized_Beta_th = ans * ax;
}

/* 分布関数: Gaussian */
double PMC::dist(AndoLab::Vector3d <double> r){

  const double Difference_of_altitude{
    std::abs( r.abs() - pR_center_magnitude )
  };

  const double Distance_in_altitude_plane{
    pR_center_magnitude * angle_between( r, pR_center )
  };

  return (
      pN0 *
      std::exp(
          - Difference_of_altitude*Difference_of_altitude
          / 2.0 / pSigma_z / pSigma_z
          - Distance_in_altitude_plane * Distance_in_altitude_plane
          / 2.0 / pSigma_r / pSigma_r )
      );
}
