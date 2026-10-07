/*
 * output.cpp
 *
 *  Created on: 2025/03/04
 *      Author: ando
 */

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>

#include "pmc.h"
#include "pmc_observation.h"
#include "pmc_simulation.h"

namespace{
std::string output_filename(const std::string &stem, const std::string &extension){
  if ( output_file_id.empty() ){
    return stem + extension;
  }
  return stem + "_" + output_file_id + extension;
}

std::string output_filename_with_index(const std::string &stem, const int index){
  if ( output_file_id.empty() ){
    return stem + "_" + std::to_string(index) + ".dat";
  }
  return stem + "_" + output_file_id + "_" + std::to_string(index) + ".dat";
}

double truncated_normal_radius_distribution(
    const double radius,
    const double mean_radius,
    const double radius_width){

  const double K0 =
      0.5 * ( 1.0 + std::erf(mean_radius / std::sqrt(2.0) / radius_width) );
  const double x = (radius - mean_radius) / radius_width;
  return std::exp(-0.5 * x * x) / std::sqrt(2.0 * M_PI) / radius_width / K0;
}

double radius_cubed_moment(
    const double mean_radius,
    const double radius_width){

  constexpr int Num_intervals { 2000 };
  const double upper_radius =
      std::max(radius_max, mean_radius + 10.0 * radius_width);
  const double step = upper_radius / static_cast <double>(Num_intervals);

  double moment = 0.0;
  for(int i = 0; i <= Num_intervals; i++){
    const double radius = step * static_cast <double>(i);
    const double weight = ( i == 0 || i == Num_intervals ) ? 1.0 : ( i % 2 == 0 ? 2.0 : 4.0 );
    moment += weight * radius * radius * radius *
        truncated_normal_radius_distribution(radius, mean_radius, radius_width);
  }

  return moment * step / 3.0;
}
}

void output_optimized_pmc_intensity(double *SolarRayIntensity, double *BackgroundIntensity,
    double **rayleigh_intensity_fine,
    std::string data_dir
    ){

  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
    std::ofstream ofs_pmc_layer(data_dir + "/" + output_filename_with_index("optimized_pmc_narrow", j_lambda));
    for(int k = 0; k < Num_PMC_Layer; k++){
      ofs_pmc_layer << (Lower_PMC + dAlt*k) * 1e-3 << " "
          << SolarRayIntensity[j_lambda]*rayleigh_intensity_fine[j_lambda][k] + BackgroundIntensity[j_lambda] << " "
          << "\n";
    }
    ofs_pmc_layer.close();
  }
}

void output_rayleigh_pmc_height(double *SolarRayIntensity, double *BackgroundIntensity,
    double **rayleigh_intensity_fine, double **rayleigh_intensity_coarse,
    std::string data_dir){

  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
    std::ofstream ofs_pmc_layer(data_dir + "/" + output_filename_with_index("pmc_rayleigh", j_lambda));
    for(int k = 0; k < Num_PMC_Layer; k++){
      ofs_pmc_layer << (Lower_PMC + dAlt*k) * 1e-3 << " "
          << SolarRayIntensity[j_lambda]*rayleigh_intensity_fine[j_lambda][k] + BackgroundIntensity[j_lambda] << " "
          << SolarRayIntensity[j_lambda]*rayleigh_intensity_coarse[j_lambda][k] + BackgroundIntensity[j_lambda] << "\n";
    }
    ofs_pmc_layer.close();
  }
}

void output_fitted(const int i_alpha,
    double **Rayleigh,
    double ***Observed_data,
    double *SolarRayIntensity, double *BackgroundIntensity,
    std::string data_dir, double &Err23_40){

  Err23_40 = 0.0;

  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
    std::ofstream ofs( data_dir + "/" + output_filename_with_index("fitting", j_lambda));
    for(int i = 0; i < N_alt; i++){
      double Sim23_40 = SolarRayIntensity[j_lambda]*Rayleigh[j_lambda][i] + BackgroundIntensity[j_lambda];
      double Obs23_40 = Observed_data[j_lambda][i_alpha][ int(Altitude_min/dAlt) + i ];
      ofs << (Altitude_min + i*dAlt) * 1e-3 << " "
          << Sim23_40 << " "
          << Obs23_40 << " "
          << "\n";
      Err23_40 += std::log(Sim23_40/Obs23_40) * std::log(Sim23_40/Obs23_40);
    }
    ofs.close();
  }

  std::ofstream ofs_err(data_dir + "/" + output_filename("err", ".dat"));
  ofs_err << Err23_40 << std::endl; /* [-] : Square Error between Fitted simulation and Observation */
  ofs_err.close();

}

/* フィッティングに用いる観測データを出力する */
void output_observed_for_fitting(const int i_alpha,
    double ***Observed_data, double **Observed_sigma, std::string data_dir){

  for(int jLambda = 0; jLambda < Num_Lambda; jLambda++){
    std::ofstream ofs_obs( data_dir + "/" + output_filename_with_index("obs", jLambda));
    for(int iAlt = 0; iAlt < Num_observed_data_altitude; iAlt++){
      ofs_obs << iAlt << " "
          << Observed_data[jLambda][i_alpha][iAlt] << " "
          << Observed_sigma
          << std::endl;
    }
    ofs_obs.close();
  }
}

void output_optimizedresult(PMC *pmc, double shift_distance, std::string data_dir,
    double latitude, double longitude, double Upper_edge){

  /*********** tangential pointにおける高さ方向を積分したPMC数(面積あたり) **********
   *
   * shift_distance を Ds とする。tangential pointから手前方向に Ds だけ
   * 回転させたところに正規分布の中心があるので、tangential pointでは
   * exp( - Ds^2 / {2σh^2} ) になる。従って、tangential pointでの、
   * 垂直方向に積分した、単位面積あたりのPMC数 Nt は
   * Nt = N0 exp( - Ds^2 / {2σh^2} ) ∫ exp( - {r-r_0}^2 / {2σz^2} ) dr
   * = N0 exp( - Ds^2 / {2σh^2} ) √2 σz ∫ exp( - t^2 ) dt
   * = N0 exp( - Ds^2 / {2σh^2} ) √(2π) σz
   */
  const double Number_of_PMC_per_area_at_tangential_point
  = pmc[0].N0() * std::exp( - shift_distance * shift_distance / 2.0 / pmc[0].sig_r() / pmc[0].sig_r() )
  * std::sqrt(2.0*M_PI) * pmc[0].sig_z();

  /*
   * Ice Water Content (IWC):
   *
   * IWC = N0 exp( - Ds^2 / {2 σh^2} ) ∫ exp( - {z - z0}^2 / { 2σz^2 } ) dz
   *       × ρ 4π/3 * ∫ r^3 f(r; r0, Δr) dr
   *      ( [m^-3] * [m] × [kg/m^3] * m^3 m^-1 * m = kg/m^2 )
   *     = N0 exp( - Ds^2 / {2 σh^2} ) √(2π) σz
   *       × ρ 4π/3 * <r^3> [kg/m^2]
   */
  const double r0 = pmc[0].r0();
  const double radius_width = pmc[0].sigma();
  const double rho = 1000.0; /* [kg/m^3] : mass density of water */
  const double IWC = pmc[0].N0()
      * std::exp( - shift_distance * shift_distance / 2.0 / pmc[0].sig_r() / pmc[0].sig_r() )
      * std::sqrt(2.0*M_PI) * pmc[0].sig_z()
      * rho * 4.0 * M_PI / 3.0 * radius_cubed_moment(r0, radius_width);

  const std::string optimized_output_file =
      output_file_id.empty() ? "output.dat" : "opt_" + output_file_id + ".dat";
  std::ofstream ofs_output( data_dir + "/" + optimized_output_file );
  ofs_output << latitude << " " << longitude << " "
      << ( pmc[0].rc().abs() - Radius_of_Earth ) * m2km << " " /* [km] : Center Altitude */
      << pmc[0].sig_z() * m2km << " " /* [km] : SD in altitude */
      << pmc[0].sig_r() * m2km << " " /* [km] : SD in horizontal direction */
      << pmc[0].r0() * 1e9 << " " /* [nm] : Particle radius */
      << pmc[0].sigma() * 1e9 << " "    /* [nm] : Width of Gaussian radius distribution */
      << pmc[0].N0() * 1e-6 << " " /* [cm^-3] : Number density at the maximum (at the center) */
      << shift_distance * m2km << " " /* [km] : horizontal shift */
      << Number_of_PMC_per_area_at_tangential_point << " " /* [m^-2] */
      << IWC * 1.0e6 * 1.0e3 * 1.0e-2 * 1.0e-2 << " " /* [kg/m^2]=>[μg/cm^2]  */
      << Upper_edge * m2km << " " /* [km] : Upper edge altitude of PMC */
      << std::endl;
  ofs_output.close();
}
