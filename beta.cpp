/*
 * 変更履歴
 * - 2026/6/17: sigma_t は標準大気で十分(組成が同一と仮定するなら、分子一つ当たりの散乱断面積は同じ)なので
 *              書き換え
 */

#include <physical_constants.h>
#include <Vector3d.h>
#include "Geocoordinate.h"
#include <Msis21.h>
#include "pmc_simulation.h"

/* 波長λ[m]のときの標準(乾燥)大気の屈折率 ー １
 *
 * from [1]
 *
 */
double ns_1(double lambda){
  lambda *= 1.0e6; /* convert from [m] to [μm]  */
  const double lambda_m2 = 1.0/lambda/lambda; /* λ^{-2} */

  /* Bodhaine Eq.(1) (Elde'n) */
  // double index
  // = 6432.8 + 2949810.0 / (146.0 - lambda_m2)
  // + 25540. / (41.0 - lambda_m2);

  /* Tomasi 2005, Eq.(17) (from Ciddor) */
  const double index
    = 5792105.0 / ( 238.0185 - lambda_m2) + 167917.0 / ( 57.362 - lambda_m2 );

  return index / 1.0e8;
}

/* King因子 from Tomasi
 * C_CO2 = 450ppmv
 */
double f(double lambda){
  lambda *= 1e6; /* λ[μm] */
  const double lambda_2 = 1.0/lambda/lambda;
  double lambda_4 = lambda_2 * lambda_2;
  constexpr double C_CO2 { 450e-6 }; /* 450 ppmv */

  return ( (0.78084 * (1.034 + 3.17e-4*lambda_2)
      + 0.20946 * (1.096 + 1.385e-3*lambda_2 + 1.448e-4*lambda_4)
            + 0.00934 + C_CO2 * 1.15 ) / (0.99964 + C_CO2)  );
}

/* depolarization factor ρ
 * King factor から求める。
 * Bucholtz1995 の Table1 と比較して差がないことを確認
 */
double rho(const double lambda){
  return 6.*( f(lambda) - 1.0 ) / ( 7.0*f(lambda) + 3.0 );
}

double gam(const double lambda){
  return rho(lambda) / ( 2.0 - rho(lambda) );
}

/*
 * 位相関数
 */
double Phase_func(const double Theta, const double lambda){
  return 3./4. * ( 1. + 3.*gam(lambda) + (1. - gam(lambda))
      * std::cos(Theta)*std::cos(Theta) )
      / ( 1.0 + 2.0 * gam(lambda) );

}

//double sigma_t(const double lambda, AndoLab::Msis21 *msis){
double sigma_t(const double lambda){

  constexpr double kB { BOLTZMANN_CONSTANT };
  constexpr double Pressure_of_standard_air { 101325. }; /* [Pa] */
  constexpr double Temperature_of_standard_air { 273.15 + 15. }; /* [K] */
  constexpr double Ns = Pressure_of_standard_air / kB / Temperature_of_standard_air; /* Number density of standard air */
  const double lambda2 = lambda*lambda;
  const double lambda4 = lambda2*lambda2;

  // const double n_1 = (msis->p() / msis->T())
  //     / (Pressure_of_standard_air / Temperature_of_standard_air)
  //     * ns_1(lambda);
  const double ns1 = ns_1(lambda);
  const double n = ns1 + 1.0;

  double sigma = 24.0 * M_PI * M_PI * M_PI
    / lambda4 / Ns / Ns
    * ns1 * ns1 * (n + 1.0) * (n + 1.0)
      / (n*n + 2.0) / (n*n + 2.0)
      * f(lambda);

  return sigma;
}

/*
 * 微分散乱断面積
 */
double sigma(const double th, const double lambda){
  return sigma(th, lambda, sigma_t(lambda));
}

double sigma(const double th, const double lambda, const double sigma_total){
  return sigma_total/4.0/M_PI * Phase_func(th, lambda);
}

/*
 * 体積散乱係数
 */
double beta(AndoLab::Vector3d <double> r, void* p){
  Parameters *param = (Parameters*)p;

  Geocoordinate gr( r );

  param->ptr_msis->calc_at(gr.altitude(), gr.latitude(), gr.longitude());
//  ofs_log << param->msis.latitude() << " "
//      << param->msis.longitude() << " "
//      << param->msis.altitude()*1e-3 << " "
//      << param->msis.N() << " "
//      << param->msis.T() << " "
//      << param->msis.p() << "\n";

  return param->sigma_total * param->ptr_msis->N();
}
