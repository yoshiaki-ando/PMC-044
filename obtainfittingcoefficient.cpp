/*
 * obtainfittingcoefficient.cpp
 *
 *  Created on: 2024/06/12
 *      Author: ando
 */
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <nlopt.hpp>

#include "pmc_simulation.h"
#include "pmc_observation.h"



class PackedParameters{
private:
public:
  double Offset;
  double *observed;
  double *simulation;
};

double SquareError(const std::vector <double> &Coef, std::vector <double> &grad, void *Params){

  /*
   * 二乗誤差を計算する
   * 観測データを Obs(i)
   * 計算データを Cal(i)とすれば、
   * E = Σ_i { log( Obs(i) ) - log( Coef * Cal(i) + Offset ) }^2 を最小化する
   *
   * ∂E/∂(Coef)
   *  = Σ_i 2 { log( Obs(i) ) - log( Coef * Cal(i) + Offset ) } { - Cal(i)/( Coef * Cal(i) + Offset ) }
   *  = -2 Σ_i { log( Obs(i) ) - log( Coef * Cal(i) + Offset ) } * Cal(i) / ( Coef * Cal(i) + Offset )
   */

  PackedParameters *p = (PackedParameters*)Params;

  double err = 0.0;
  double grad_err = 0.0;
  for(int i = 0; i < N_alt; i++){
    double logdiff = std::log(p->observed[i]) - std::log(Coef[0] * p->simulation[i] + p->Offset);
    err += logdiff * logdiff;
    grad_err += logdiff * p->simulation[i] / ( Coef[0] * p->simulation[i] + p->Offset );
  }

  if ( !grad.empty() ){
    grad[0] = -2.0 * grad_err;
  }

  number_of_iteration++;
  if ( number_of_iteration > 10000 ){
    throw nlopt::forced_stop();
  }


  return err / N_alt;
}

void ObtainFittingCoefficient(
    double **Rayleigh,          /* シミュレーションしたデータ */
    const int i_alpha,          /* 観測データのうち、使用する緯度・経度のインデックス */
    double ***Observed_data,    /* 観測データ */
    double *SolarRayIntensity,  /* 係数 = 太陽光の強さ */
    double *BackgroundIntensity /* オフセット = 背景光の明るさ */
    ){

  nlopt::opt opt( nlopt::LN_NELDERMEAD, 1 );

  for(int jLambda = 0; jLambda < Num_Lambda; jLambda++){ /* 波長毎に評価は全く別 */

    /* 背景強度を、高度(Idx_background)km以上から平均で計算 */
    double AvgBgIntensity = 0.0;
    for(int iAlt = Idx_background; iAlt < Num_observed_data_altitude; iAlt++){
      AvgBgIntensity += Observed_data[jLambda][i_alpha][iAlt];
    }
    BackgroundIntensity[jLambda] = AvgBgIntensity / ( (Num_observed_data_altitude-1) - Idx_background + 1 );
    /* ↑ Note: Num_observed_data_altitude-1 が最大のインデックス */

    PackedParameters p;
    p.Offset = BackgroundIntensity[jLambda];
    const int idx_min = std::round( Altitude_min / dAlt );
//    std::cout << "Index min. = " << idx_min << " " << Altitude_min / dAlt << std::endl;

    p.observed = Observed_data[jLambda][i_alpha] + idx_min;
    p.simulation = Rayleigh[jLambda];


    opt.set_min_objective(SquareError, (void*)(&p) );
    opt.set_xtol_rel(1.0e-2);
    std::vector <double> x(1, 1.0e4);
    double minf;

    try {
      number_of_iteration = 0;
      nlopt::result result = opt.optimize(x, minf);
    } catch (std::exception &e){
      std::cout << "Optimization of fitting coefficient (" << process_id << ") failed. : " << e.what() << std::endl;
    }

//    std::cout << jLambda << ", Optimized: " << x[0] << std::endl;
    SolarRayIntensity[jLambda] = x[0];

  }

}
