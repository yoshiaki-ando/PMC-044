/*
 * msis_result.h
 *
 * MSISを用いてレイリー散乱の結果を記憶するためのクラス
 *
 *  Created on: 2024/10/16
 *      Author: ando
 */

#ifndef MSIS_RESULT_H_
#define MSIS_RESULT_H_

#include <Vector3d.h>
#include "Geocoordinate.h"

constexpr int REGION1 { 0 };
constexpr int REGION2 { 1 };
constexpr int REGION3 { 2 };
constexpr int REGION4 { 3 };
constexpr int REGION5 { 4 };

enum class CALC_MSIS { NO, YES, RAYLEIGH_INTVL, PMC_INTVL };
/* MSISを計算しない（保存したものを用いる、区間分割は自動で割り振る）、
 * MSISを計算する（区間分割は自動で割り振る）
 * 強制的にRayleigh散乱のみの区間分割で計算するか、
 * 強制的にPMC用の区間分割で計算するか */

class Msis_Result{
  /* Region毎の変数が必要 */
private:
  bool allocate_regions_done { false };
public:
  int Number_of_Regions;
  double **OpticalDepth_LOS;
  double **OpticalDepth_ToScatteringPoint;

  AndoLab::Vector3d <double> *length; /* 積分区間全長 */
  int *Num; /* 分割点数 */
  AndoLab::Vector3d <double> *dr; /* 分割した区間のベクトル */
  double *dr_abs;

  double **total_beta;
  double **N;

  AndoLab::Vector3d <double> **r_p;
  AndoLab::Vector3d <double> **r_i;

  bool **is_illuminated;

  void allocate_regions(void){
    if ( allocate_regions_done ){
      return;
    }

    OpticalDepth_LOS               = new double* [Number_of_Regions];
    OpticalDepth_ToScatteringPoint = new double* [Number_of_Regions];

    length = new AndoLab::Vector3d<double> [Number_of_Regions];
    Num = new int [Number_of_Regions];
    dr = new AndoLab::Vector3d<double> [Number_of_Regions];
    dr_abs = new double [Number_of_Regions];

    total_beta = new double* [Number_of_Regions];
    N = new double* [Number_of_Regions];

    r_p = new AndoLab::Vector3d<double>* [Number_of_Regions];
    r_i = new AndoLab::Vector3d<double>* [Number_of_Regions];

    is_illuminated = new bool* [Number_of_Regions];

    allocate_regions_done = true;
  };

  void allocate_dr(const int Region_num){
    OpticalDepth_LOS[Region_num]               = new double [Num[Region_num]+1];
    OpticalDepth_ToScatteringPoint[Region_num] = new double [Num[Region_num]];
    total_beta[Region_num] = new double[Num[Region_num]];
    N[Region_num] = new double [Num[Region_num]];

    r_p[Region_num] = new AndoLab::Vector3d <double> [ Num[Region_num] ];
    r_i[Region_num] = new AndoLab::Vector3d <double> [ Num[Region_num] ];

    is_illuminated[Region_num] = new bool [ Num[Region_num] ];
  }
};

#endif /* MSIS_RESULT_H_ */
