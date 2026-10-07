/*
 * レイリー散乱の計算
 */
#include <iostream>
#include <fstream>
#include <array>
#include <string>

#include "Fixed_parameter.h"
#include "Across_pts.h"
#include "pmc_simulation.h"

extern std::ofstream ofs_log;

constexpr int j_pmc { 0 }; /* 将来的な複数PMC対応のための予備 */

void calculate_intensity(
    Date Day_of_Year,
    const double lambda,
    const double Lower_Altitude,
    const double Step_Altitude,
    const int Number_of_Altitude,
    const double alpha,
    AndoLab::Msis21 &msis,
    PMC &pmc,
    double *intensity,             /* 返り値。各高度の光強度 */
    CALC_MSIS calc_msis,
    Msis_Result **msis_result
){

  std::string function_name( "[calculate_intensity] " );

  /* ガウス積分に渡すクラス化したパラメタ
   * 本当はもっと上で定義できる */
  Parameters param;
  param.ptr_msis = &msis;  /* MSISE */
  param.wavelength = lambda;
  param.sigma_total = sigma_t(lambda);

  AndoLab::Vector3d <double> Solar_d = r_s(Day_of_Year); /* 太陽の方向 */

  for(int iAlt = 0; iAlt < Number_of_Altitude; iAlt++){

    /* 計算結果を保存するための記憶確保 */


//    std::cout << iAlt << " / " << Number_of_Altitude << " : " << std::endl;
//    ofs_log << iAlt << " / " << Number_of_Altitude << " : " << std::endl;
    intensity[iAlt] = 0.0;
    double alt = Lower_Altitude + iAlt*Step_Altitude;

    /* 高度alt, 北極からの角度αのtangential point */
    AndoLab::Vector3d <double> r_tangential = tangential_point(alt, alpha);
    /* 衛星からの視線と、大気圏上界との交点 */
    std::array<AndoLab::Vector3d <double>, 2> Pts_atmos = Across_point_atmosphere(r_tangential);


    AndoLab::Vector3d <double> Pts_upper_pmc[2], Pts_lower_pmc[2];
    /* PMC領域上面との交点 */
    if ( alt < Upper_PMC ){
      Across_pts PtsAcrossUpperPMC(Pts_atmos[0], Pts_atmos[1], Upper_PMC);
      PtsAcrossUpperPMC.copy(Pts_upper_pmc);
    }
    /* PMC領域下面との交点 */
    if ( alt < Lower_PMC ){
      Across_pts PtsAcrossLowerPMC(Pts_atmos[0], Pts_atmos[1], Lower_PMC);
      PtsAcrossLowerPMC.copy(Pts_lower_pmc);
    }

    /* 散乱角
     * Pts_atmos[0] - Pts_atmos[1] は、衛星の方向
     * -Solar_d は太陽光の入射方向
     */
    double th =
        angle_between( Pts_atmos[0] - Pts_atmos[1], -1.0*Solar_d);
    param.sigma_differential = sigma(th, lambda, param.sigma_total);

    /* Mie散乱のための前処理(2) */
    if ( pmc.is_this_real() ){
        pmc.calc_normalized_beta_th(th);
    }


    /* 以下で、各高度の散乱光強度を計算。高速化の必要あり */
    if ( alt < Lower_PMC ){

      if ( calc_msis != CALC_MSIS::NO ){
        const int Number_of_Regions { 5 };

        for(int m_intvl = 0; m_intvl < 2; m_intvl++){
          msis_result[iAlt][m_intvl].Number_of_Regions = Number_of_Regions;
          msis_result[iAlt][m_intvl].allocate_regions();
        }
      }

      double delta2r = 0.0;
      /* (1) 大気圏上部からPMC層上部までの計算
       * Pts_atmos[0]     -> Pts_upper_pmc[0] …散乱点はレイリー */
      intensity[iAlt] +=
          intensity_integral(Pts_atmos[0], Pts_upper_pmc[0],
              Day_of_Year, Solar_d, lambda,
              th, msis, delta2r, pmc, param, HIGHER_THAN_PML_LAYER, Integral_Interval,
              calc_msis, REGION1, msis_result[iAlt]);

      /* (2) PMC層
       * Pts_upper_pmc[0] -> Pts_lower_pmc[0] …散乱点はPMC層内 */
//      std::cout << "region2" << std::endl;
      if ( (pmc.is_this_real() &&  This_position_is_near_PMC( pmc, Pts_upper_pmc[0], Pts_lower_pmc[0] ))
          || (calc_msis == CALC_MSIS::PMC_INTVL) ){
        /* PMCがあるときは細かな刻みで積分する */
        intensity[iAlt] +=
            intensity_integral(Pts_upper_pmc[0], Pts_lower_pmc[0],
                Day_of_Year, Solar_d, lambda,
                th, msis, delta2r, pmc, param, INSIDE_PML_LAYER, PMC_Integral_Interval,
                calc_msis, REGION2, msis_result[iAlt]);

      } else {
        intensity[iAlt] +=
            intensity_integral(Pts_upper_pmc[0], Pts_lower_pmc[0],
                Day_of_Year, Solar_d, lambda,
                th, msis, delta2r, pmc, param, INSIDE_PML_LAYER, Integral_Interval,
                calc_msis, REGION2, msis_result[iAlt]);
      }

      /* (3) PMC層下部から、PMC層下部まで
       * Pts_lower_pmc[0] -> Pts_lower_pmc[1] …散乱点はレイリー */
//      std::cout << "region3" << std::endl;
      intensity[iAlt] +=
          intensity_integral(Pts_lower_pmc[0], Pts_lower_pmc[1],
              Day_of_Year, Solar_d, lambda,
              th, msis, delta2r, pmc, param, LOWER_THAN_PML_LAYER, Integral_Interval,
              calc_msis, REGION3, msis_result[iAlt]);

      /* (4) PMC層
       * Pts_lower_pmc[1] -> Pts_upper_pmc[1] …散乱点はPMC層内 */
//      std::cout << "region4" << std::endl;
      if ( (pmc.is_this_real() &&  This_position_is_near_PMC( pmc, Pts_lower_pmc[1], Pts_upper_pmc[1] ))
          || (calc_msis == CALC_MSIS::PMC_INTVL) ){
        /* PMCがあるときは細かな刻みで積分する */
        intensity[iAlt] +=
            intensity_integral(Pts_lower_pmc[1], Pts_upper_pmc[1],
                Day_of_Year, Solar_d, lambda,
                th, msis, delta2r, pmc, param, INSIDE_PML_LAYER, PMC_Integral_Interval,
                calc_msis, REGION4, msis_result[iAlt]);

      } else {
        intensity[iAlt] +=
            intensity_integral(Pts_lower_pmc[1], Pts_upper_pmc[1],
                Day_of_Year, Solar_d, lambda,
                th, msis, delta2r, pmc, param, INSIDE_PML_LAYER, Integral_Interval,
                calc_msis, REGION4, msis_result[iAlt]);

      }

      /* (5) PMC層上部から大気圏上部まで
       * Pts_upper_pmc[1] -> Pts_atmos[1]     …散乱点はレイリー */
//      std::cout << "region5" << std::endl;
      intensity[iAlt] +=
          intensity_integral(Pts_upper_pmc[1], Pts_atmos[1],
              Day_of_Year, Solar_d, lambda,
              th, msis, delta2r, pmc, param, HIGHER_THAN_PML_LAYER, Integral_Interval,
              calc_msis, REGION5, msis_result[iAlt]);


    } else if ( alt < Upper_PMC ){ /* PMC内 */

      if ( calc_msis != CALC_MSIS::NO ){
        const int Number_of_Regions { 3 };

        for(int m_intvl = 0; m_intvl < 2; m_intvl++){
          msis_result[iAlt][m_intvl].Number_of_Regions = Number_of_Regions;
          msis_result[iAlt][m_intvl].allocate_regions();
        }
      }


      double delta2r = 0.0;
      /* (1) 大気圏上部からPMC層上部までの計算
       * Pts_atmos[0]     -> Pts_upper_pmc[0] …散乱点はレイリー */
      intensity[iAlt] +=
          intensity_integral(Pts_atmos[0], Pts_upper_pmc[0],
              Day_of_Year, Solar_d, lambda,
              th, msis, delta2r, pmc, param, HIGHER_THAN_PML_LAYER, Integral_Interval,
              calc_msis, REGION1, msis_result[iAlt]);

      /* (2) PMC層
       * Pts_upper_pmc[0] -> Pts_upper_pmc[1] …散乱点はPMC層内 */
      if ( (pmc.is_this_real() &&  This_position_is_near_PMC( pmc, Pts_upper_pmc[0], Pts_upper_pmc[1] ))
          || (calc_msis == CALC_MSIS::PMC_INTVL) ){

        /* PMCがあるときは細かな刻みで積分する */
        intensity[iAlt] +=
            intensity_integral(Pts_upper_pmc[0], Pts_upper_pmc[1],
                Day_of_Year, Solar_d, lambda,
                th, msis, delta2r, pmc, param, INSIDE_PML_LAYER, PMC_Integral_Interval,
                calc_msis, REGION2, msis_result[iAlt]);
      } else {

        intensity[iAlt] +=
            intensity_integral(Pts_upper_pmc[0], Pts_upper_pmc[1],
                Day_of_Year, Solar_d, lambda,
                th, msis, delta2r, pmc, param, INSIDE_PML_LAYER, Integral_Interval,
                calc_msis, REGION2, msis_result[iAlt]);
      }

      /* (3) PMC層上部から大気圏上部まで
       * Pts_upper_pmc[1] -> Pts_atmos[1]   …散乱点はレイリー   */
      intensity[iAlt] +=
          intensity_integral(Pts_upper_pmc[1], Pts_atmos[1],
              Day_of_Year, Solar_d, lambda,
              th, msis, delta2r, pmc, param, HIGHER_THAN_PML_LAYER, Integral_Interval,
              calc_msis, REGION3, msis_result[iAlt]);

    } else {

      if ( calc_msis != CALC_MSIS::NO  ){
        const int Number_of_Regions { 1 };

        for(int m_intvl = 0; m_intvl < 2; m_intvl++){
          msis_result[iAlt][m_intvl].Number_of_Regions = Number_of_Regions;
          msis_result[iAlt][m_intvl].allocate_regions();
        }
      }

      double delta2r = 0.0;
      /* Pts_atmos[0] -> Pts_atmos[1]   …散乱点はレイリー   */
      intensity[iAlt] =
          intensity_integral(Pts_atmos[0], Pts_atmos[1],
              Day_of_Year, Solar_d, lambda,
              th, msis, delta2r, pmc, param, HIGHER_THAN_PML_LAYER, Integral_Interval,
              calc_msis, REGION1, msis_result[iAlt]);

    }

//    std::cout << std::endl;
  }

}
