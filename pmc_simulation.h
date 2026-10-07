/*
 * pmc_simulation.h
 *
 *  Created on: 2024/05/27
 *      Author: ando
 */

#ifndef PMC_SIMULATION_H_
#define PMC_SIMULATION_H_

#include <array>
#include <string>

#include <Msis21.h>
#include <Vector3d.h>

#include "Fixed_parameter.h"
#include "Date.h"
#include "pmc.h"
#include "msis_result.h"


/* 衛星視線上で大気圏上界面間を、分けて積分する区間 */
constexpr double Integral_Interval { 200.0e3 };
constexpr double PMC_Integral_Interval { 1.0e3 };


/*********************
 * 設定するパラメタ *
 *********************/

/* フィッティングパラメタを計算する高度 */
constexpr double Altitude_min { 20.e3 };
constexpr double Altitude_max { 43.e3 };
constexpr int N_alt { 24 }; /* 43 - 20 + 1 = 24 */
//constexpr double Altitude_min { 0.e3 };
//constexpr double Altitude_max { 100.e3 };
//constexpr int N_alt { 101 }; /* 43 - 20 + 1 = 24 */
constexpr double dAlt { (Altitude_max - Altitude_min) / (N_alt-1) };

constexpr int Maximum_convergence { 800 }; /* PMC最適化での反復回数上限 */

/*******************************
 * 設定するパラメタ (ここまで) *
 *******************************/

enum class LOGGING { NO, COLLECTED, SUCCESSIVE };
extern LOGGING logging_method;
extern std::string logging_str;

class Parameters{
private:
public:
  double wavelength;
  double sigma_total;
  double sigma_differential;
  AndoLab::Msis21 *ptr_msis;
};

/* 収束しないときの反復回数 */
extern int number_of_iteration;

/* ジョブID */
extern std::string process_id;

/* 出力ファイル名に付ける識別子 YYYYMMDD_hhmm_latlon-index */
extern std::string output_file_id;

/********************** 関数プロトタイプ **********************/

/* 引数から緯度、経度、日付を取得 */
void get_arg(const int argc, char **argv, double &latitude, double &longitude,
    Date &date,
    AndoLab::Msis21 &msis, std::string &data_dir,
    std::string &process_id);

/* 緯度・経度から、計算に使う座標(ひまわり方向を +x方向)へ変換 */
void convert_coordinate(
    const double latitude, /* [deg] */
    const double longitude, /* [deg] */
    AndoLab::Vector3d <double> &r0 /* [m] converted coordinate at altitude 0km */
    );

void calculate_intensity(
    Date Day_of_Year,
    const double Lambda,
    const double Lower_Altitude,
    const double Step_Altitude,
    const int Number_of_Altitude,
    const double alpha,
    AndoLab::Msis21 &msis,
    PMC &pmc,
    double *intensity,             /* 返り値。各高度の光強度 */
    CALC_MSIS calculate_msis_or_not,
    Msis_Result **msis_result /* Altitude X PMC/Rayleigh_interval */
    );

/* ひまわりから見て、北極から角度αで、高度Altitudeの球面とのtangential point
 * αが正なら東経140°以下側、負なら東経140°以上 */
AndoLab::Vector3d <double> tangential_point(double Altitude_in_m, double Angle_from_North_Pole_in_rad);

/* ひまわりから r の点を見た時、大気圏上空(高度100km)との交点2点を計算 */
std::array<AndoLab::Vector3d <double>, 2> Across_point_atmosphere(AndoLab::Vector3d <double> r);

/* 点r から d方向にいったとき、高度 H で交差する点 */
AndoLab::Vector3d <double> cross_point_at_altitude(
    AndoLab::Vector3d <double> r,
    AndoLab::Vector3d <double> d, /* 単位ベクトルとすること */
    const double H /* [m] 高度 */
    );

/* 太陽を示すベクトル */
AndoLab::Vector3d <double> r_s(Date);

/*
 * r1 → r2 の経路中、pmc の分布と近くなる点が存在するか
 */
bool This_position_is_near_PMC(PMC &pmc, AndoLab::Vector3d <double> &r1, AndoLab::Vector3d <double> &r2 );

double intensity_integral(
    AndoLab::Vector3d <double> &r_from,
    AndoLab::Vector3d <double> &r_to,
    Date &date,
    const AndoLab::Vector3d <double> &Solar_d,
    const double Wavelength,
    const double th, /* 散乱角 */
    AndoLab::Msis21 &msis,
    double &delta,
    PMC &pmc,
    Parameters &param,
    const int Whether_inside_or_outside_PMC_region,
    const double integral_interval,
    CALC_MSIS calculate_msis_or_not,
    const int Region_number,
    Msis_Result *msis_result
    );

double beta(AndoLab::Vector3d <double> r, void* p);
double sigma_t(const double lambda);
double sigma(const double th, const double lambda);
double sigma(const double th, const double lambda, const double sigma_total);

inline double lat2theta(double latitude){
  return (90 - latitude)*Deg2rad;
}
inline double phi(double theta){
  return acos( Radius_of_Earth / Rgeo / sin(theta) );
}

void ObtainFittingPMC(
    Date date,
    PMC *pmc,
    AndoLab::Msis21 &msis,
    const double alpha,
    const int i_alpha,
    double ***Observed_data,
    double **rayleigh_intensity,
    const double *SolarRayIntensity,
    const double *BackgroundIntensity,
    double &shift_distance,
    Msis_Result ***msis_result,
    std::string data_dir
    );

std::complex <double> interpolated_m(double wavelength_in_m);

void output_optimizedresult(PMC *pmc, double shift_distance, std::string data_dir,
    double latitude, double longitude, double Upper_edge);
void output_observed_for_fitting(const int i_alpha,
    double ***Observed_data, double **Observed_sigma, std::string data_dir);
void output_fitted(const int i_alpha,
    double **Rayleigh,
    double ***Observed_data,
    double *SolarRayIntensity, double *BackgroundIntensity,
    std::string data_dir, double &Err23_40);
void output_rayleigh_pmc_height(double *SolarRayIntensity, double *BackgroundIntensity,
    double **rayleigh_intensity_fine, double **rayleigh_intensity_coarse,
    std::string data_dir);
void output_optimized_pmc_intensity(double *SolarRayIntensity, double *BackgroundIntensity,
    double **rayleigh_intensity_fine,
    std::string data_dir );


void solar_direction(
 const double lat,
 const double lon,
 const int doy,
 const double hours,
 double &A,
 double &th_h,
 double *ZenithAngle,
 double *Azimuth
 );
#endif /* PMC_SIMULATION_H_ */
