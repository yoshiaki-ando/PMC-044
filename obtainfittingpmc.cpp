/*
 * obtainfittingpmc.cpp
 *
 *  Created on: 2024/06/18
 *      Author: ando
 */
#include <iostream>
#include <fstream>
#include <vector>
//#include <chrono>
#include <string>
#include <array>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <nlopt.hpp>

#include <Msis21.h>
#include <memory_allocate.h>

#include "pmc_simulation.h"
#include "msis_result.h"


enum class OptParam : int {
  IdxAltitude = 0,
      IdxSigmaZ = 1,
      IdxSigmaH = 2,
      IdxR0 = 3,
      IdxSigma = 4,
      IdxShift = 5,
      IdxNc = 6
};

constexpr int NUM_OPTIMIZED_PARAMETER { 7 };
constexpr double MAX_PMC_CENTER_DENSITY { 1000.0e6 };

std::vector <double> normalize_parameters(
    const std::vector <double> &physical_param,
    const std::vector <double> &LowerBound,
    const std::vector <double> &UpperBound){

  std::vector <double> normalized_param(NUM_OPTIMIZED_PARAMETER);
  for(int i = 0; i < NUM_OPTIMIZED_PARAMETER; i++){
    normalized_param[i] =
        ( physical_param[i] - LowerBound[i] ) / ( UpperBound[i] - LowerBound[i] );
  }

  return normalized_param;
}

std::vector <double> denormalize_parameters(
    const std::vector <double> &normalized_param,
    const std::vector <double> &LowerBound,
    const std::vector <double> &UpperBound){

  std::vector <double> physical_param(NUM_OPTIMIZED_PARAMETER);
  for(int i = 0; i < NUM_OPTIMIZED_PARAMETER; i++){
    physical_param[i] =
        LowerBound[i] + normalized_param[i] * ( UpperBound[i] - LowerBound[i] );
  }

  return physical_param;
}

double column_density_from_center_density(
    const double center_density,
    const double sigma_z,
    const double sigma_h,
    const double shift){

  return center_density
      * std::exp( - shift * shift / 2.0 / sigma_h / sigma_h )
      * std::sqrt(2.0 * M_PI) * sigma_z;
}

double center_density_from_column_density(
    const double column_density,
    const double sigma_z,
    const double sigma_h,
    const double shift){

  const double column_per_center_density =
      column_density_from_center_density(1.0, sigma_z, sigma_h, shift);
  if ( column_per_center_density <= 0.0 ){
    return std::numeric_limits <double>::infinity();
  }

  return column_density / column_per_center_density;
}

double center_density_from_optimized_param(const std::vector <double> &Optimized_param){

  return center_density_from_column_density(
      Optimized_param[static_cast <int>(OptParam::IdxNc)],
      Optimized_param[static_cast <int>(OptParam::IdxSigmaZ)],
      Optimized_param[static_cast <int>(OptParam::IdxSigmaH)],
      Optimized_param[static_cast <int>(OptParam::IdxShift)]);
}

AndoLab::Vector3d <double> rotate_alpha(AndoLab::Vector3d <double> r, const double alpha);
AndoLab::Vector3d <double> rotate_theta(AndoLab::Vector3d <double> r, const double theta);

AndoLab::Vector3d <double> set_center_pmc(
    const double pmc_altitude, const double horizontal_shift, const double alpha){

  const double R = Radius_of_Earth + pmc_altitude;
  double cos_b = R / Rgeo;
  double sin_b = std::sqrt( 1.0 - cos_b*cos_b );
  AndoLab::Vector3d <double> r(R*cos_b, 0.0, R*sin_b );
  return rotate_alpha( rotate_theta(r, horizontal_shift/R ), alpha );
}

extern std::ofstream ofs_log;

class PMC_FittingParameter{
private:
public:
  int i_alpha;
  double alpha;
  PMC *pmc;
  const Date *ptr_date;
  AndoLab::Msis21 *ptr_msis;
  const double *SolarRayIntensity;
  const double *BackgroundIntensity;
  double ***Observed_data;
  double **rayleigh_intensity;
  Msis_Result ***msis_result;
  const std::vector <double> *LowerBound;
  const std::vector <double> *UpperBound;
};

//constexpr double N0_base { 20.0e6 };

/* 同定に用いる評価関数。この誤差を最小化する */
double ErrPMC(const std::vector <double> &Normalized_param, std::vector <double> &grad, void *Params){

//  std::chrono::system_clock::time_point start_time
//  = std::chrono::system_clock::now();

  PMC_FittingParameter *p = (PMC_FittingParameter*)Params;
  const std::vector <double> Optimized_param =
      denormalize_parameters(Normalized_param, *(p->LowerBound), *(p->UpperBound));
  const double center_density = center_density_from_optimized_param(Optimized_param);
  if ( !std::isfinite(center_density)
      || center_density < 0.0
      || center_density > MAX_PMC_CENTER_DENSITY ){
    number_of_iteration++;
    return 1.0e100;
  }

  double **intensity = AndoLab::allocate_memory2d(Num_Lambda, Num_PMC_Layer, 0.0);

  /* PMCの中心 */
  AndoLab::Vector3d <double> rc =
      set_center_pmc(Optimized_param[static_cast <int>(OptParam::IdxAltitude)],
          Optimized_param[static_cast <int>(OptParam::IdxShift)], p->alpha); /* PMCの中心 */

  double sq_err { 0.0 }; /* 二乗誤差 */

  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
    std::complex<double> m_ice = interpolated_m( Lambda[j_lambda] );
    p->pmc[j_lambda].set( Lambda[j_lambda], m_ice, rc,
                          Optimized_param[static_cast <int>(OptParam::IdxSigmaZ)],
                          Optimized_param[static_cast <int>(OptParam::IdxSigmaH)],
                          Optimized_param[static_cast <int>(OptParam::IdxR0)],
                          Optimized_param[static_cast <int>(OptParam::IdxSigma)],
                          center_density );
    /* PMCの高度のみを求める */
    calculate_intensity( *( p->ptr_date ), Lambda[j_lambda],
        Lower_PMC, dAlt, Num_PMC_Layer,
        p->alpha, *(p->ptr_msis), p->pmc[j_lambda], intensity[j_lambda],
        CALC_MSIS::NO, p->msis_result[j_lambda] );

//    std::ofstream ofs("data/searching_" + std::to_string(number_of_iteration) + "_"
//        + std::to_string(j_lambda) + ".dat");
//    for(int k = 0; k < Num_PMC_Layer; k++){
//      ofs << (Lower_PMC + dAlt*k) * 1e-3 << " "
//          << p->SolarRayIntensity[j_lambda]*intensity[j_lambda][k] +
//          p->BackgroundIntensity[j_lambda] << "\n";
//    }
//    ofs.close();
  }

  /* 数密度係数を求める
   *
   */
//  double nume { 0.0 };
//  double deno { 0.0 };
//  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
//    double nume_lambda { 0.0 };
//    double deno_lambda { 0.0 };
//    for(int idx_alt = 0; idx_alt < Num_PMC_Layer; idx_alt++){
//      double I_PMC = intensity[j_lambda][idx_alt] - p->rayleigh_intensity[j_lambda][idx_alt];
//      nume_lambda += I_PMC * ( p->Observed_data[j_lambda][p->i_alpha][IdxLowerAlt + idx_alt] -
//          p->SolarRayIntensity[j_lambda] * p->rayleigh_intensity[j_lambda][idx_alt] -
//          p->BackgroundIntensity[j_lambda] );
//      deno_lambda += I_PMC * I_PMC;
//    }
//    nume += p->SolarRayIntensity[j_lambda] * nume_lambda;
//    deno += p->SolarRayIntensity[j_lambda] * p->SolarRayIntensity[j_lambda] * deno_lambda;
//  }
//  double density_coefficient { nume / deno };
//  std::cout << density_coefficient << std::endl;

//    std::ofstream ofs("data/diff_" + std::to_string(j_lambda) + ".dat");
  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
    /* 対数値の二乗誤差の計算 */
    for(int idx_alt = 0; idx_alt < Num_PMC_Layer; idx_alt++){
//      double I_PMC = intensity[j_lambda][idx_alt] - p->rayleigh_intensity[j_lambda][idx_alt];
//      double calc_result = p->SolarRayIntensity[j_lambda] *
//          ( p->rayleigh_intensity[j_lambda][idx_alt] + density_coefficient * I_PMC ) +
//          p->BackgroundIntensity[j_lambda];
      double calc_result =
          p->SolarRayIntensity[j_lambda] * intensity[j_lambda][idx_alt] +
          p->BackgroundIntensity[j_lambda];
      const int IdxLowerAlt { int( std::round( Lower_PMC / dAlt ) ) };
      double log_diffs =
          std::log( p->Observed_data[j_lambda][p->i_alpha][IdxLowerAlt + idx_alt] / calc_result );
      sq_err += log_diffs * log_diffs;
    }
    //    ofs.close();
  }
  double SqErr = std::sqrt( sq_err );

  AndoLab::deallocate_memory2d( intensity );

//  p->pmc[0].N0( density_coefficient * N0_base ); /* 数密度係数で、N0を設定 */

  /* 最適化している様子を出力 */
  std::string str_log = "(" + process_id + ") "
      + std::to_string( number_of_iteration ) + " "
      + std::to_string( SqErr ) + " " /* Error */
      + std::to_string( ( p->pmc[0].rc().abs() - Radius_of_Earth ) * m2km ) + " " /* Center Altitude */
      + std::to_string( p->pmc[0].sig_z() * m2km ) + " " /* SD in altitude */
      + std::to_string( p->pmc[0].sig_r() * m2km ) + " " /* SD in horizontal */
      + std::to_string( p->pmc[0].r0() * 1e9 ) + " " /* Particle radius */
      + std::to_string( p->pmc[0].sigma() * 1e9 ) + " "    /* Width of Gaussian radius distribution */
      + std::to_string( Optimized_param[static_cast <int>(OptParam::IdxNc)] ) + " " /* Column density at tangential point */
      + std::to_string( p->pmc[0].N0() * 1e-6 ) + " " /* Number density at the maximum (at the center) */
      + std::to_string(  Optimized_param[static_cast <int>(OptParam::IdxShift)] * m2km ) + " "; /* horizontal shift */

  if ( logging_method == LOGGING::COLLECTED ){
    logging_str += str_log + "\n";
  } else if ( logging_method == LOGGING::SUCCESSIVE ){
    ofs_log << str_log << std::endl;
  }

//  std::chrono::system_clock::time_point end_time
//  = std::chrono::system_clock::now();
//  double elapsed
//     = std::chrono::duration_cast <std::chrono::milliseconds>
//     (end_time - start_time).count();
//  std::cout << elapsed * 1e-3 << " [sec]" << std::endl;
//  ofs_log << elapsed * 1e-3 << " [sec]" << std::endl;

  number_of_iteration++;
//  if ( number_of_iteration > Maximum_convergence ){
//    throw nlopt::forced_stop();
//  }
  return SqErr;
}

double evaluate_pmc_error_without_log(
    const std::vector <double> &physical_param,
    PMC_FittingParameter &Param,
    const std::vector <double> &LowerBound,
    const std::vector <double> &UpperBound){

  std::vector <double> normalized_param =
      normalize_parameters(physical_param, LowerBound, UpperBound);
  std::vector <double> grad(NUM_OPTIMIZED_PARAMETER);

  const LOGGING previous_logging_method = logging_method;
  const int previous_number_of_iteration = number_of_iteration;

  logging_method = LOGGING::NO;
  const double error = ErrPMC(normalized_param, grad, (void*)(&Param));
  logging_method = previous_logging_method;
  number_of_iteration = previous_number_of_iteration;

  return error;
}

void output_one_parameter_sweep(
    PMC_FittingParameter &Param,
    const std::vector <double> &Optimized_param,
    const std::vector <double> &LowerBound,
    const std::vector <double> &UpperBound,
    std::string data_dir){

  constexpr int Num_sweep_points { 41 };
  const std::array <std::string, NUM_OPTIMIZED_PARAMETER> ParamName {
    "altitude", "sigma_z", "sigma_h", "r0", "radius_width", "shift", "Nc"
  };
  const std::array <std::string, NUM_OPTIMIZED_PARAMETER> UnitName {
    "km", "km", "km", "nm", "nm", "km", "m^-2"
  };
  const std::array <double, NUM_OPTIMIZED_PARAMETER> OutputScale {
    m2km, m2km, m2km, 1.0e9, 1.0e9, m2km, 1.0
  };

  const double base_error =
      evaluate_pmc_error_without_log(Optimized_param, Param, LowerBound, UpperBound);

  for(int idx_param = 1; idx_param < NUM_OPTIMIZED_PARAMETER; idx_param++){
    const std::string filename =
        data_dir + "/sweep_" + output_file_id + "_" + ParamName[idx_param] + ".dat";
    std::ofstream ofs(filename);
    ofs << std::setprecision(16);
    ofs << "# " << ParamName[idx_param] << " [" << UnitName[idx_param] << "] "
        << "error delta_error relative_delta_error\n";
    ofs << "# base_error " << base_error << "\n";

    for(int i = 0; i < Num_sweep_points; i++){
      const double ratio =
          static_cast <double>(i) / static_cast <double>(Num_sweep_points - 1);
      std::vector <double> candidate_param = Optimized_param;
      candidate_param[idx_param] =
          LowerBound[idx_param] + ratio * ( UpperBound[idx_param] - LowerBound[idx_param] );

      const double error =
          evaluate_pmc_error_without_log(candidate_param, Param, LowerBound, UpperBound);
      ofs << candidate_param[idx_param] * OutputScale[idx_param] << " "
          << error << " "
          << error - base_error << " "
          << ( error - base_error ) / base_error << "\n";
    }
  }
}

class PMC_OptimizationResult{
private:
public:
  std::string stage;
  int index;
  int evaluations;
  double error;
  std::vector <double> physical_param;
};

void append_optimizer_log_comment(const std::string &comment){

  if ( logging_method == LOGGING::COLLECTED ){
    logging_str += "# " + comment + "\n";
  } else if ( logging_method == LOGGING::SUCCESSIVE ){
    ofs_log << "# " << comment << std::endl;
  }
}

void output_optimization_summary_header(std::ofstream &ofs){

  ofs << "# stage index evaluations error "
      << "CenterAlt[km] Sig_z[km] Sig_r[km] R0[nm] RadiusWidth[nm] "
      << "Shift[km] Nc[m^-2] N0[x10^6_m^-3]\n";
}

void output_optimization_summary(
    std::ofstream &ofs,
    const PMC_OptimizationResult &result){

  ofs << result.stage << " "
      << result.index << " "
      << result.evaluations << " "
      << result.error << " "
      << result.physical_param[static_cast <int>(OptParam::IdxAltitude)] * m2km << " "
      << result.physical_param[static_cast <int>(OptParam::IdxSigmaZ)] * m2km << " "
      << result.physical_param[static_cast <int>(OptParam::IdxSigmaH)] * m2km << " "
      << result.physical_param[static_cast <int>(OptParam::IdxR0)] * 1.0e9 << " "
      << result.physical_param[static_cast <int>(OptParam::IdxSigma)] * 1.0e9 << " "
      << result.physical_param[static_cast <int>(OptParam::IdxShift)] * m2km << " "
      << result.physical_param[static_cast <int>(OptParam::IdxNc)] << " "
      << center_density_from_optimized_param(result.physical_param) * 1.0e-6 << "\n";
}

int getenv_int_or_default(const char *name, const int default_value){

  const char *value = std::getenv(name);
  if ( value == nullptr ){
    return default_value;
  }

  const int int_value = std::atoi(value);
  if ( int_value <= 0 ){
    return default_value;
  }

  return int_value;
}

PMC_OptimizationResult run_local_pmc_optimizer(
    PMC_FittingParameter &Param,
    const std::vector <double> &initial_physical_param,
    const std::vector <double> &LowerBound,
    const std::vector <double> &UpperBound,
    const std::vector <double> &NormalizedLowerBound,
    const std::vector <double> &NormalizedUpperBound,
    const std::vector <double> &InitialStep,
    const nlopt::algorithm optimizer_algorithm,
    const int maxeval,
    const std::string &stage,
    const int index){

  const int initial_iteration = number_of_iteration;
  std::vector <double> normalized_param =
      normalize_parameters(initial_physical_param, LowerBound, UpperBound);

  append_optimizer_log_comment(stage + " index=" + std::to_string(index)
      + " maxeval=" + std::to_string(maxeval));

  nlopt::opt opt(optimizer_algorithm, NUM_OPTIMIZED_PARAMETER);
  opt.set_min_objective(ErrPMC, (void*)(&Param));
  opt.set_xtol_rel(1.0e-2);
  opt.set_initial_step(InitialStep);
  opt.set_lower_bounds(NormalizedLowerBound);
  opt.set_upper_bounds(NormalizedUpperBound);
  opt.set_maxeval(maxeval);

  double minf = std::numeric_limits <double>::infinity();
  try {
    opt.optimize(normalized_param, minf);
  } catch (std::exception &e){
    std::cout << "Optimization of PMC parameters (" << process_id << ", "
        << stage << " " << index << ") failed. : " << e.what() << std::endl;
    append_optimizer_log_comment("Optimization failed in " + stage
        + " index=" + std::to_string(index) + " : " + e.what());
  }

  PMC_OptimizationResult result;
  result.stage = stage;
  result.index = index;
  result.evaluations = number_of_iteration - initial_iteration;
  result.physical_param = denormalize_parameters(normalized_param, LowerBound, UpperBound);
  result.error = evaluate_pmc_error_without_log(
      result.physical_param, Param, LowerBound, UpperBound);

  return result;
}


void ObtainFittingPMC(
    /*
     ********* PMCの最適化 *********
     * NLOpt を用いて PMCパラメタを同定する
     */
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
    std::string data_dir){


  /*
   * Optimized_param[0] : [m] 高度 Lower_PMC - Upper_PMC
   * Optimized_param[1] : [m] σ_z 0.1km - 5.0km
   * Optimized_param[2] : [m] σ_h (horizontal) 10km - 2000km
   * Optimized_param[3] : [m] r0 粒径分布の平均半径 10nm - 180nm
   * Optimized_param[4] : [m] Δr 粒径分布の幅 2nm - 30nm
   * Optimized_param[5] : [m] r_s  +y方向に対して右ねじにPMC中心位置をこの距離だけ回転(r_s = (R0 + z_pmc)*θ )
   * Optimized_param[6] : [m^-2] Nc 接線点における柱粒子数密度
   */
  std::vector <double> Optimized_param(NUM_OPTIMIZED_PARAMETER);
  Optimized_param[static_cast <int>(OptParam::IdxAltitude)] = 82.5057e3;
  Optimized_param[static_cast <int>(OptParam::IdxSigmaZ)] = 1.10469e3;
  Optimized_param[static_cast <int>(OptParam::IdxSigmaH)] = 750.0e3;
  Optimized_param[static_cast <int>(OptParam::IdxR0)] = 66.5536e-9;
  Optimized_param[static_cast <int>(OptParam::IdxSigma)] = 11.3e-9;
  Optimized_param[static_cast <int>(OptParam::IdxShift)] = 150e3;
  Optimized_param[static_cast <int>(OptParam::IdxNc)] =
      column_density_from_center_density(22.3356e6,
          Optimized_param[static_cast <int>(OptParam::IdxSigmaZ)],
          Optimized_param[static_cast <int>(OptParam::IdxSigmaH)],
          Optimized_param[static_cast <int>(OptParam::IdxShift)]);

  const char *initial_opt_file = std::getenv("PMC_INITIAL_OPT_FILE");
  if ( initial_opt_file != nullptr ){
    std::ifstream ifs_initial(initial_opt_file);
    double latitude_dummy, longitude_dummy;
    double altitude_km, sigma_z_km, sigma_h_km, r0_nm, radius_width_nm;
    double N0_cm3, shift_km;
    if ( ifs_initial >> latitude_dummy >> longitude_dummy
        >> altitude_km >> sigma_z_km >> sigma_h_km >> r0_nm
        >> radius_width_nm >> N0_cm3 >> shift_km ){
      Optimized_param[static_cast <int>(OptParam::IdxAltitude)] = altitude_km * 1.0e3;
      Optimized_param[static_cast <int>(OptParam::IdxSigmaZ)] = sigma_z_km * 1.0e3;
      Optimized_param[static_cast <int>(OptParam::IdxSigmaH)] = sigma_h_km * 1.0e3;
      Optimized_param[static_cast <int>(OptParam::IdxR0)] = r0_nm * 1.0e-9;
      Optimized_param[static_cast <int>(OptParam::IdxSigma)] = radius_width_nm * 1.0e-9;
      Optimized_param[static_cast <int>(OptParam::IdxShift)] = shift_km * 1.0e3;
      double column_density_m2;
      if ( ifs_initial >> column_density_m2 ){
        Optimized_param[static_cast <int>(OptParam::IdxNc)] = column_density_m2;
      } else {
        Optimized_param[static_cast <int>(OptParam::IdxNc)] =
            column_density_from_center_density(N0_cm3 * 1.0e6,
                Optimized_param[static_cast <int>(OptParam::IdxSigmaZ)],
                Optimized_param[static_cast <int>(OptParam::IdxSigmaH)],
                Optimized_param[static_cast <int>(OptParam::IdxShift)]);
      }
      std::cout << "PMC initial opt file : " << initial_opt_file << std::endl;
    } else {
      std::cout << "Could not read PMC initial opt file : " << initial_opt_file << std::endl;
    }
  }

  /* 最大値・最小値の設定 */
  std::vector <double> LowerBound { Lower_PMC, 0.1e3, 10.0e3,    10.0e-9, 2.0e-9, 0.0, 0.0 };
  std::vector <double> UpperBound { Upper_PMC, 5.0e3, 2000.0e3, 180.0e-9, 30.0e-9, 1500.0e3, 1.0e12 };
  std::vector <double> NormalizedLowerBound(NUM_OPTIMIZED_PARAMETER, 0.0);
  std::vector <double> NormalizedUpperBound(NUM_OPTIMIZED_PARAMETER, 1.0);
  std::vector <double> InitialStep(NUM_OPTIMIZED_PARAMETER, 0.1);

  /* わたすパラメタの設定 */
  PMC_FittingParameter Param;
  Param.alpha = alpha;
  Param.pmc = pmc;
  Param.ptr_date = &date;
  Param.ptr_msis = &msis;
  Param.SolarRayIntensity = SolarRayIntensity;
  Param.BackgroundIntensity = BackgroundIntensity;
  Param.Observed_data = Observed_data;
  Param.i_alpha = i_alpha;
  Param.rayleigh_intensity = rayleigh_intensity;
  Param.msis_result = msis_result;
  Param.LowerBound = &LowerBound;
  Param.UpperBound = &UpperBound;

//  ErrPMC( Optimized_param, grad, (void*)(&Param) ); /* お試し計算 */

  std::ofstream ofs_optimization_summary(
      data_dir + "/optimization_summary_" + output_file_id + ".dat");
  ofs_optimization_summary << std::setprecision(16);
  output_optimization_summary_header(ofs_optimization_summary);

  number_of_iteration = 0;

  /* NLOpt による最適化 */
  bool use_staged_optimizer { true };
  nlopt::algorithm optimizer_algorithm = nlopt::LN_BOBYQA;
  std::string optimizer_name { "staged LN_BOBYQA" };
  const char *optimizer_env = std::getenv("PMC_OPTIMIZER");
  if ( optimizer_env != nullptr ){
    const std::string requested_optimizer { optimizer_env };
    if ( requested_optimizer == "NELDERMEAD" || requested_optimizer == "LN_NELDERMEAD" ){
      optimizer_algorithm = nlopt::LN_NELDERMEAD;
      optimizer_name = "LN_NELDERMEAD";
      use_staged_optimizer = false;
    } else if ( requested_optimizer == "BOBYQA" || requested_optimizer == "LN_BOBYQA" ){
      optimizer_algorithm = nlopt::LN_BOBYQA;
      optimizer_name = "LN_BOBYQA";
      use_staged_optimizer = false;
    } else if ( requested_optimizer == "STAGED" || requested_optimizer == "MULTISTART" ){
      optimizer_algorithm = nlopt::LN_BOBYQA;
      optimizer_name = "staged LN_BOBYQA";
      use_staged_optimizer = true;
    }
  }
  if ( initial_opt_file != nullptr && optimizer_env == nullptr ){
    optimizer_algorithm = nlopt::LN_BOBYQA;
    optimizer_name = "LN_BOBYQA";
    use_staged_optimizer = false;
  }

  std::cout << "PMC optimizer : " << optimizer_name << std::endl;
  append_optimizer_log_comment("PMC optimizer: " + optimizer_name);

  if ( use_staged_optimizer ){
    const int Stage1Maxeval = getenv_int_or_default("PMC_STAGE1_MAXEVAL", 40);
    const int Stage2Maxeval = getenv_int_or_default("PMC_STAGE2_MAXEVAL", 500);
    const int Stage3Maxeval = getenv_int_or_default("PMC_STAGE3_MAXEVAL", 1000);
    constexpr int NumStage2Starts { 3 };

    const std::array <double, 2> Stage1SigmaH_km { 300.0, 1300.0 };
    const std::array <double, 2> Stage1R0_nm { 30.0, 60.0 };
    const std::array <double, 2> Stage1Shift_km { 100.0, 700.0 };
    const std::array <double, 2> Stage1N0_cm3 { 20.0, 100.0 };

    // const std::array <double, 4> Stage1SigmaH_km { 300, 800, 1300, 1800 };
    // const std::array <double, 3> Stage1R0_nm { 40, 70, 100 };
    // const std::array <double, 3> Stage1Shift_km { 100.0, 300.0, 1000.0 };
    // const std::array <double, 3> Stage1N0_cm3 { 20, 50, 150 };

    // const std::array <double, 4> Stage1SigmaH_km { 100, 500, 1900, 1200 };
    // const std::array <double, 4> Stage1R0_nm { 20, 30, 55, 80 };
    // const std::array <double, 3> Stage1Shift_km { 200.0, 700.0, 1400.0 };
    // /* Stage1N0_cm3 is used only to build the initial Nc. The optimized parameter is Nc. */
    // const std::array <double, 4> Stage1N0_cm3 { 10, 40, 70, 100 };

    std::vector <PMC_OptimizationResult> stage1_results;
    int stage1_index = 0;
    for(double sigma_h_km : Stage1SigmaH_km){
      for(double r0_nm : Stage1R0_nm){
        for(double shift_km : Stage1Shift_km){
          for(double N0_cm3 : Stage1N0_cm3){
            std::vector <double> initial_param(NUM_OPTIMIZED_PARAMETER);
            initial_param[static_cast <int>(OptParam::IdxAltitude)] = 82.0e3;
            initial_param[static_cast <int>(OptParam::IdxSigmaZ)] = 1.0e3;
            initial_param[static_cast <int>(OptParam::IdxSigmaH)] = sigma_h_km * 1.0e3;
            initial_param[static_cast <int>(OptParam::IdxR0)] = r0_nm * 1.0e-9;
            initial_param[static_cast <int>(OptParam::IdxSigma)] = 15.0e-9;
            initial_param[static_cast <int>(OptParam::IdxShift)] = shift_km * 1.0e3;
            initial_param[static_cast <int>(OptParam::IdxNc)] =
                column_density_from_center_density(N0_cm3 * 1.0e6,
                    initial_param[static_cast <int>(OptParam::IdxSigmaZ)],
                    initial_param[static_cast <int>(OptParam::IdxSigmaH)],
                    initial_param[static_cast <int>(OptParam::IdxShift)]);

            PMC_OptimizationResult result =
                run_local_pmc_optimizer(Param, initial_param,
                    LowerBound, UpperBound, NormalizedLowerBound, NormalizedUpperBound,
                    InitialStep, nlopt::LN_BOBYQA, Stage1Maxeval, "stage1", stage1_index);
            stage1_results.push_back(result);
            output_optimization_summary(ofs_optimization_summary, result);
            stage1_index++;
          }
        }
      }
    }

    std::sort(stage1_results.begin(), stage1_results.end(),
        [](const PMC_OptimizationResult &a, const PMC_OptimizationResult &b){
          return a.error < b.error;
        });

    for(int i = 0; i < NumStage2Starts; i++){
      PMC_OptimizationResult top_result = stage1_results[i];
      top_result.stage = "stage1_top";
      top_result.index = i;
      output_optimization_summary(ofs_optimization_summary, top_result);
    }

    std::vector <PMC_OptimizationResult> stage2_results;
    for(int i = 0; i < NumStage2Starts; i++){
      PMC_OptimizationResult result =
          run_local_pmc_optimizer(Param, stage1_results[i].physical_param,
              LowerBound, UpperBound, NormalizedLowerBound, NormalizedUpperBound,
              InitialStep, nlopt::LN_BOBYQA, Stage2Maxeval, "stage2", i);
      stage2_results.push_back(result);
      output_optimization_summary(ofs_optimization_summary, result);
    }

    std::sort(stage2_results.begin(), stage2_results.end(),
        [](const PMC_OptimizationResult &a, const PMC_OptimizationResult &b){
          return a.error < b.error;
        });

    PMC_OptimizationResult final_result =
        run_local_pmc_optimizer(Param, stage2_results[0].physical_param,
            LowerBound, UpperBound, NormalizedLowerBound, NormalizedUpperBound,
            InitialStep, nlopt::LN_BOBYQA, Stage3Maxeval, "stage3", 0);
    output_optimization_summary(ofs_optimization_summary, final_result);
    Optimized_param = final_result.physical_param;
  } else {
    PMC_OptimizationResult result =
        run_local_pmc_optimizer(Param, Optimized_param,
            LowerBound, UpperBound, NormalizedLowerBound, NormalizedUpperBound,
            InitialStep, optimizer_algorithm, 1000, "single", 0);
    output_optimization_summary(ofs_optimization_summary, result);
    Optimized_param = result.physical_param;
  }

  const double optimized_center_density = center_density_from_optimized_param(Optimized_param);

  /* pmcクラスに残っているパラメタ */
  std::string str_log = "# pmc class\n"
      + std::to_string(number_of_iteration) + " "
      + std::to_string( ( pmc[0].rc().abs() - Radius_of_Earth ) * m2km ) + " " /* Center Altitude [km] */
      + std::to_string(pmc[0].sig_z() * m2km) + " " /* SD in altitude [km] */
      + std::to_string(pmc[0].sig_r() * m2km) + " " /* SD in horizontal [km] */
      + std::to_string( pmc[0].r0() * 1e9 ) + " " /* Particle radius [nm] */
      + std::to_string(pmc[0].sigma() * 1e9 ) + " " /* Width of Gaussian radius distribution [nm] */
      + std::to_string(Optimized_param[static_cast <int>(OptParam::IdxNc)] ) + " " /* Column density at tangential point [m^-2] */
      + std::to_string(pmc[0].N0() * 1e-6) + " " /* Number density at the maximum (at the center) [x10^6 m^-3] */
      + "\n";
//  ofs_log << "# pmc class\n"
//      << number_of_iteration  << " "
//      << ( pmc[0].rc().abs() - Radius_of_Earth ) * m2km << " " /* Center Altitude [km] */
//      << pmc[0].sig_z() * m2km << " " /* SD in altitude [km] */
//      << pmc[0].sig_r() * m2km << " " /* SD in horizontal [km] */
//      << pmc[0].r0() * 1e9 << " " /* Particle radius [nm] */
//      << pmc[0].sigma() * 1e9 << " "    /* Width of Gaussian radius distribution [nm] */
//      << Optimized_param[static_cast <int>(OptParam::IdxNc)] << " " /* Column density at tangential point [m^-2] */
//      << pmc[0].N0() * 1e-6 << " " /* Number density at the maximum (at the center) [x10^6 m^-3] */
//      << std::endl;

  /* 最適化されたパラメタ */
  str_log += "# opt param\n"
      + std::to_string(Optimized_param[static_cast <int>(OptParam::IdxAltitude)] * m2km) + " "
      + std::to_string(Optimized_param[static_cast <int>(OptParam::IdxSigmaZ)] * m2km) + " "
      + std::to_string(Optimized_param[static_cast <int>(OptParam::IdxSigmaH)] * m2km) + " "
      + std::to_string(Optimized_param[static_cast <int>(OptParam::IdxR0)] * 1e9) + " "
      + std::to_string(Optimized_param[static_cast <int>(OptParam::IdxSigma)] * 1e9 ) + " "
      + std::to_string(Optimized_param[static_cast <int>(OptParam::IdxShift)] * m2km) + " "
      + std::to_string(Optimized_param[static_cast <int>(OptParam::IdxNc)] ) + " "
      + std::to_string(optimized_center_density * 1e-6) + " ";

//  ofs_log << "# opt param\n"
//      << Optimized_param[static_cast <int>(OptParam::IdxAltitude)] * m2km << " "
//      << Optimized_param[static_cast <int>(OptParam::IdxSigmaZ)] * m2km << " "
//      << Optimized_param[static_cast <int>(OptParam::IdxSigmaH)] * m2km << " "
//      << Optimized_param[static_cast <int>(OptParam::IdxR0)] * 1e9 << " "
//      << Optimized_param[static_cast <int>(OptParam::IdxSigma)] * 1e9 << " "
//      << Optimized_param[static_cast <int>(OptParam::IdxShift)] * m2km << " "
//      << Optimized_param[static_cast <int>(OptParam::IdxNc)] << " "
//      << optimized_center_density * 1e-6 << " "
//      << std::endl;

  if ( logging_method == LOGGING::COLLECTED ){
    logging_str += str_log + "\n";
  } else if ( logging_method == LOGGING::SUCCESSIVE ){
    ofs_log << str_log << std::endl;
  }

  if ( std::getenv("PMC_SWEEP") != nullptr ){
    output_one_parameter_sweep(Param, Optimized_param, LowerBound, UpperBound, data_dir);
  }

  /* 最適化された値をPMCに設定する */
//  AndoLab::Vector3d <double> rc = tangential_point( Optimized_param[0], alpha );
  AndoLab::Vector3d <double> rc =
      set_center_pmc(Optimized_param[static_cast <int>(OptParam::IdxAltitude)],
          Optimized_param[static_cast <int>(OptParam::IdxShift)], alpha);

//  const double N0 = pmc[0].N0();
  double *intensity = new double [Num_PMC_Layer]; /* デバッグ追加 */

//  PMC dummy_pmc;

  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){

    std::complex<double> m_ice = interpolated_m( Lambda[j_lambda] );
    pmc[j_lambda].set( Lambda[j_lambda], m_ice, rc,
        Optimized_param[static_cast <int>(OptParam::IdxSigmaZ)], /* σ_z */
        Optimized_param[static_cast <int>(OptParam::IdxSigmaH)], /* σ_h */
        Optimized_param[static_cast <int>(OptParam::IdxR0)], /* r0 */
        Optimized_param[static_cast <int>(OptParam::IdxSigma)], /* Δr */
        optimized_center_density); /* N0 */


    calculate_intensity(date, Lambda[j_lambda],
        Lower_PMC, dAlt, Num_PMC_Layer,
//        alpha, msis, dummy_pmc, intensity,
        alpha, msis, pmc[j_lambda], intensity,
        CALC_MSIS::NO, msis_result[j_lambda]);  /* デバッグ追加 */

//    calculate_intensity(date, Lambda[j_lambda],
//        Lower_PMC, dAlt, Num_PMC_Layer,
//        alpha, msis, dummy_pmc, intensity2,
//        alpha, msis, pmc[j_lambda], intensity2,
//        CALC_MSIS::YES, msis_result[j_lambda]);  /* デバッグ追加 */
//        CALC_MSIS::YES, msis_optimized[j_lambda]);  /* デバッグ追加 */
    std::ofstream ofs( data_dir + "/optimized_pmc_insubroutine_" + output_file_id + "_"
        + std::to_string(j_lambda) + ".dat"); /* デバッグ追加 */
    for(int i = 0; i < Num_PMC_Layer; i++){ /* デバッグ追加 */
      ofs << Lower_PMC*1e-3 + i << " "
          << SolarRayIntensity[j_lambda]*intensity[i] + BackgroundIntensity[j_lambda] << "\n";
    } /* デバッグ追加 */
    ofs.close(); /* デバッグ追加 */

  }
//  delete [] msis_optimized[0];
//  delete [] msis_optimized[1];
//  delete [] msis_optimized[2];
//  delete [] msis_optimized;

  shift_distance = Optimized_param[static_cast <int>(OptParam::IdxShift)];


  delete [] intensity; /* デバッグ追加 */

}
