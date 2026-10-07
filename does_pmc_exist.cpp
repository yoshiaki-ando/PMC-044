/*
 * does_pmc_exist.cpp
 *
 *  Created on: 2024/12/03
 *      Author: ando
 */
#include <iostream>
#include <fstream>
#include "Fixed_parameter.h"
#include "pmc_simulation.h"
#include "pmc_observation.h"

bool does_pmc_exist(
    const int i_alpha,
    double ***Observed_Average,
    double **Observed_SD,
    double **Simulated_Rayleigh,    /* 観測データ */
    double *SolarRayIntensity,  /* 係数 = 太陽光の強さ */
    double *BackgroundIntensity, /* オフセット = 背景光の明るさ */
    double &Upper_edge /* PMCがあると判定された最高高度 */
    ){

  bool exists = false;
  Upper_edge = -1.0;

  constexpr int jLambda { 0 }; /* 青の光のみ用いて判定する */

  /* 背景光σ */
  double AvgBgSigma = 0.0;
  for(int iAlt = Idx_background; iAlt < Num_observed_data_altitude; iAlt++){
    AvgBgSigma += Observed_SD[i_alpha][iAlt];
  }
  AvgBgSigma /= ( Num_observed_data_altitude - Idx_background );

  int iAlt = idx_UpperPMC;
  /* レイリー散乱 + 背景光平均σ < 観測光平均 - 0.5*観測光σ
   * が成立するならPMCが存在する */
  while ( (iAlt >= idx_LowerPMC) && !exists ){
//  while ( (iAlt <= idx_UpperPMC) ){
    double Simulated_Intensity =
        SolarRayIntensity[jLambda] * Simulated_Rayleigh[jLambda][iAlt-idx_LowerPMC] + BackgroundIntensity[jLambda];
//    ofs << iAlt << " "
//        << Simulated_Intensity << " "
//        << AvgBgSigma << " "
//        << Observed_Average[jLambda][i_alpha][iAlt] <<  " "
//        << Observed_SD[i_alpha][iAlt] << " "
//        << (Simulated_Intensity + AvgBgSigma < Observed_Average[jLambda][i_alpha][iAlt] - 0.5 * Observed_SD[i_alpha][iAlt])
//        << std::endl;
    exists = exists ||
        (Simulated_Intensity + AvgBgSigma < Observed_Average[jLambda][i_alpha][iAlt] - 0.5 * Observed_SD[i_alpha][iAlt]);
    iAlt--;
  }
//  ofs.close();

  if ( exists ){
    Upper_edge = (iAlt + 1) * DAltitude_observed;
  }

  return exists;
}


