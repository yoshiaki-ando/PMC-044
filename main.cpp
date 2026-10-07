/*
 * PMC-043
 * 計算内容：一つの緯度・経度位置におけるPMCパラメタの同定を行う
 *           PMC-041は多数の再計算があるので、計算結果をメモリに入れて高速化を図る
 *
 * 実行方法：
 *
 *
 * 出力ファイル(ディレクトリは実行ディレクトリ下の data/ )：
 *  1. rayleigh_X.dat : (Xは0-2で各波長)、フィッティングする係数を求めるために計算した
 *    レイリー散乱のデータ
 *  2. obs_X.dat : フィッティングする観測データ
 *  3. optimized_pmc_X.dat : 同定されたPMCパラメタを用いて、レイリー＋PMCによる散乱を計算
 *  4. log.txt : 最適化している途中経過。誤差、
 *
 * main.cpp 緯度 経度 YYMMDD(or YYYYMMDD) hhmm LatLonIndex
 *
 *  Created on: 2024/05/27
 *      Author: ando
 */
#include <iostream>
#include <fstream>
#include <filesystem>
#include <complex>
#include <chrono>
#include <string>

#include <Msis21.h>
#include <memory_allocate.h>
#include <Mie_scattering.h>

#include "Fixed_parameter.h"
#include "Date.h"
#include "Geocoordinate.h"
#include "pmc.h"
#include "pmc_simulation.h"
#include "pmc_observation.h"
#include "msis_result.h"
#include "write_log.h"

std::ofstream ofs_log; /* ログ出力 */
LOGGING logging_method { LOGGING::COLLECTED }; /* ログ出力の方法 */
std::string logging_str { "" };

int number_of_iteration;

std::ofstream ofs_param; /* 計算したパラメタの出力 */

std::string process_id;
std::string output_file_id;

std::string obs_data_dir("/home/ando/98-Local/PMC_data/convert_data");

int main(int argc, char **argv){

  wrap_msisinit_(); /* MSIS初期化 */

  /************************************************************
   * 初期化
   ************************************************************/

  Date date; /* (自作の日付クラス)解析をする日付、UT */
  double latitude, longitude; /* 解析をする(観測点の)緯度・経度 */

  AndoLab::Msis21 msis;
  std::string data_dir;

  /* 引数から、観測点の緯度・経度、(msisクラスへの)日付・時刻を取得 */
  get_arg(argc, argv, latitude, longitude, date, msis, data_dir, process_id);

  /* data/YYYY/YYYYMM/YYYYMMDD に出力する */
  std::filesystem::create_directories( data_dir );

  ofs_param.open( data_dir + "/param_" + output_file_id + ".txt", std::ios::out);

  write_log_open(data_dir); /* ログ出力・オープン */

  /* 緯度・経度から、計算に使う座標(ひまわり方向を +x方向)への変換
   *
   * 緯度・経度・高度 0m の位置を r0 とする。地表面のtangential point
   */
  AndoLab::Vector3d <double> r0;
  convert_coordinate(latitude, longitude, r0);
  Geocoordinate Geo_r0( r0 ); /* ひまわり座標・地理座標を扱う変換 */
  std::cout << process_id << " : "
      << latitude << ", " << longitude << " ==> "
      << Geo_r0.latitude() << ", " << Geo_r0.longitude() << std::endl;

  double alpha = Geo_r0.alpha(); /* ひまわりから見た角度 */

  /************************************************************
   * 初期化（ここまで）
   ************************************************************/

  /************************************************************
   * (TODO)太陽天頂角が90°以上なら終了
   ************************************************************/

  /************************************************************
   * 20-43kmのレイリー散乱を用いてフィッティングを行う
   ************************************************************/

  /* レイリー散乱のみ計算 */
  double **Rayleigh = AndoLab::allocate_memory2d(Num_Lambda, N_alt, 0.0);
  PMC pmc[Num_Lambda]; /* 未設定のダミーPMCで、レイリー散乱部分を計算 */
  Msis_Result ***dummy_msis_result = new Msis_Result** [Num_Lambda]; /* 積分した結果を保存(使わない) */

  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
    dummy_msis_result[j_lambda] = new Msis_Result* [N_alt];
    for(int k_alt = 0; k_alt < N_alt; k_alt++){
      dummy_msis_result[j_lambda][k_alt] = new Msis_Result [2]; /* PMC/Rayleigh interval */
    }

    calculate_intensity(date, Lambda[j_lambda],
        Altitude_min, dAlt, N_alt,
        alpha, msis, pmc[j_lambda], Rayleigh[j_lambda],
        CALC_MSIS::RAYLEIGH_INTVL, dummy_msis_result[j_lambda]);
  }


  /* 観測データの読み込み */
  /* 観測データ(平均・分散)の取得 */
  double ***Observed_data
  = AndoLab::allocate_memory3d(Num_Lambda, Num_observed_data_latitude, Num_observed_data_altitude, 0.0);
  /* 分散は波長0 (青)だけ */
  double **Observed_sigma
  = AndoLab::allocate_memory2d(Num_observed_data_latitude, Num_observed_data_altitude, 0.0);
  double **Obsrvd_LatLon = AndoLab::allocate_memory2d(2, Num_observed_data_latitude, 0.0);
  get_observation_data(Observed_data, Observed_sigma, Obsrvd_LatLon, argv);

  /* 観測データのうち、この計算で求める緯度・経度のデータのインデックスを取得する */
  const int i_alpha = select_latlon( Obsrvd_LatLon, latitude, longitude );


  /* フィッティングに用いる観測データを出力する */
  output_observed_for_fitting(i_alpha, Observed_data, Observed_sigma, data_dir);

  /* フィッティング係数の導出 */
  double *SolarRayIntensity = new double [Num_Lambda];   /* 係数 = 太陽光の強さ */
  double *BackgroundIntensity = new double [Num_Lambda]; /* オフセット = 背景光の明るさ */
  ObtainFittingCoefficient(Rayleigh, i_alpha, Observed_data, SolarRayIntensity, BackgroundIntensity);

  /* フィッティング結果の出力 */
  /*** TODO: フィッティングが出来てないものはチェックする必要がある ***/
  double Err23_40;
  output_fitted(i_alpha, Rayleigh, Observed_data, SolarRayIntensity, BackgroundIntensity,
      data_dir, Err23_40);

  /* PMC高度におけるレイリー散乱(レイリー散乱のみ求める積分区間で)を求めておく
   * 計算結果を記憶する
   */
  Msis_Result ***msis_result = new Msis_Result** [Num_Lambda];
  double **rayleigh_intensity_coarse = AndoLab::allocate_memory2d(Num_Lambda, Num_PMC_Layer, 0.0);

  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){

    msis_result[j_lambda] = new Msis_Result* [Num_PMC_Layer];
    for(int k_alt = 0; k_alt < Num_PMC_Layer; k_alt++){
      msis_result[j_lambda][k_alt] = new Msis_Result [2];
    }

    calculate_intensity(date, Lambda[j_lambda],
        Lower_PMC, dAlt, Num_PMC_Layer,
        alpha, msis, pmc[j_lambda], rayleigh_intensity_coarse[j_lambda],
        CALC_MSIS::RAYLEIGH_INTVL, msis_result[j_lambda]);
  }

  /* PMCが存在するかの判断 by Kawaura and Tsuda method */
  double Upper_edge;
  bool pmc_exists = does_pmc_exist(i_alpha, Observed_data, Observed_sigma,
      rayleigh_intensity_coarse, SolarRayIntensity, BackgroundIntensity, Upper_edge);
  std::cout << "PMC exists : " << pmc_exists << std::endl;
  if ( !pmc_exists ){ /* PMCが存在しない */
    exit(0);
  }

  /* PMC高度におけるレイリー散乱(PMCの区間で)を求めておく
   * 計算結果を記憶する
   */
  double **rayleigh_intensity_fine = AndoLab::allocate_memory2d(Num_Lambda, Num_PMC_Layer, 0.0);
  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){

    calculate_intensity(date, Lambda[j_lambda],
        Lower_PMC, dAlt, Num_PMC_Layer,
        alpha, msis, pmc[j_lambda], rayleigh_intensity_fine[j_lambda],
        CALC_MSIS::PMC_INTVL, msis_result[j_lambda]);

  }
  /* フィッティングした PMC高度(78-88) を出力 */
  output_rayleigh_pmc_height(SolarRayIntensity, BackgroundIntensity,
      rayleigh_intensity_fine, rayleigh_intensity_coarse, data_dir);

  /* PMCパラメタの最適化をする */
  double shift_distance;
  ObtainFittingPMC(date, pmc, msis, alpha, i_alpha,
      Observed_data, rayleigh_intensity_fine, SolarRayIntensity, BackgroundIntensity,
      shift_distance, msis_result, data_dir);

//  /* 最適化されたPMCパラメタを用いて、 0-99kmを計算 */
//  constexpr int Num_all_allaltitude { 100 };
//
//  double *intensity = new double [Num_all_allaltitude];
//  Msis_Result ***msis_optimized = new Msis_Result** [Num_Lambda];
//
////  output_log = true;
//
//  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
//    msis_optimized[j_lambda] = new Msis_Result* [ Num_all_allaltitude ];
//    for(int k_alt = 0; k_alt < Num_all_allaltitude; k_alt++){
//      msis_optimized[j_lambda][k_alt] = new Msis_Result [2];
//    }
//
//    calculate_intensity(date, Lambda[j_lambda],
//        0.0e3, dAlt, Num_all_allaltitude,
//        alpha, msis, pmc[j_lambda], intensity,
//        CALC_MSIS::YES, msis_optimized[j_lambda]);
//
//    std::ofstream ofs( data_dir + "/optimized_pmc_" + std::to_string(j_lambda) + ".dat");
//    for(int i = 0; i < 100; i++){
//      ofs << i << " " << SolarRayIntensity[j_lambda]*intensity[i] + BackgroundIntensity[j_lambda] << "\n";
//    }
//    ofs.close();
//    ofs_log << "Used [" << j_lambda << "] : " << SolarRayIntensity[j_lambda] << ", " << BackgroundIntensity[j_lambda] << std::endl;
//  }
//  delete [] intensity;
  /* 最適化されたPMCパラメタを用いて、 78-88kmを計算 */
  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){

    calculate_intensity(date, Lambda[j_lambda],
        Lower_PMC, dAlt, Num_PMC_Layer,
        alpha, msis, pmc[j_lambda], rayleigh_intensity_fine[j_lambda],
        CALC_MSIS::NO, msis_result[j_lambda]);
  }
  output_optimized_pmc_intensity( SolarRayIntensity, BackgroundIntensity,
      rayleigh_intensity_fine, data_dir );

  /* 最適化パラメタの出力 */
  output_optimizedresult(pmc, shift_distance, data_dir,
      latitude, longitude, Upper_edge);

  AndoLab::deallocate_memory2d(Rayleigh);
  AndoLab::deallocate_memory3d(Observed_data);
  AndoLab::deallocate_memory2d(Observed_sigma);
  AndoLab::deallocate_memory2d(Obsrvd_LatLon);
  AndoLab::deallocate_memory2d(rayleigh_intensity_fine);
  AndoLab::deallocate_memory2d(rayleigh_intensity_coarse);
  delete [] SolarRayIntensity;
  delete [] BackgroundIntensity;

  ofs_param.close();

  write_log_close(data_dir);

  return 0;
}
