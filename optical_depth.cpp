#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <stdio.h>

#include <Vector3d.h>
#include <vector_gauss_quadrature.h>
#include <Msis21.h>

#include "Across_pts.h"
#include "Date.h"
#include "Geocoordinate.h"
#include "Mie_scattering.h"
#include "pmc.h"
#include "Region.h"
#include "pmc_simulation.h"

extern std::ofstream ofs_log;


//Region search_pmc_region(
//    AndoLab::Vector3d <double> &r2,
//    AndoLab::Vector3d <double> &r1
//    );

double mie_beta(AndoLab::Vector3d <double> r, void *parameters){

  PMC *mp = (PMC*)parameters;
//  std::cout << r.x()*m2km << ", " << r.y()*m2km << ", " << r.z()*m2km << std::endl;
//  std::cout << mp->dist(r) << std::endl;
//  std::cout << mp->bn() << std::endl;

  return mp->dist(r) * mp->bn();
}

double Mie_opt_depth(
    AndoLab::Vector3d <double> &r2, AndoLab::Vector3d <double> &r1,
    PMC &pmc){

  AndoLab::Vector3d <double> d = r2 - r1;
  const double distance = d.abs();
  const int Num = int(distance / PMC_Integral_Interval) + 1;
  AndoLab::Vector3d <double> dr = d/double(Num);

  double ans = 0.0;
  for(int i = 0; i < Num; i++){
    ans += vector_gauss_quadrature(mie_beta, (void*)(&pmc), r2+double(i)*dr, r2+double(i+1)*dr);
  }

  return ans;
}

double Optical_depth(
    AndoLab::Vector3d <double> &r2, AndoLab::Vector3d <double> &r1,
    PMC &pmc, Parameters &param, const bool pmc_region){

  /* r2 -> r1 への光学的深さ */


  /* レイリー散乱による光学的深さ
   * (ガウス積分で全区間を積分してしまう→問題ない) */
//  const double distance { (r2-r1).abs() };
//  const int N_inte { int( distance / Integral_Interval ) + 1 };
//  AndoLab::Vector3d <double> dr = (r2-r1) / double(N_inte);
//  double delta { 0.0 };
//  for(int i = 0; i < N_inte; i++){
//    delta += gauss_integration_vector(beta, (void*)(&param),
//        r1 + double(i+1)*dr, r1 + double(i)*dr, gauss_z, gauss_w);
//  }
  double delta = vector_gauss_quadrature(beta, (void*)(&param), r2, r1);


  /* Mie散乱の寄与 */
  if ( pmc_region && pmc.is_this_real() ){
    /*
     * PMCの分布から遠い場合は計算しない
     * 「遠い」の判断は
     * 条件A: r1 と r2 の両方が、PMCの中心rc から(PMCの分布の標準偏差)×3以上の距離にあり、かつ
     * 条件B: (r1 - rc)・(r2 - rc) > 0 →同じ方向
     * 条件A and 条件B なら PMCの光学的深さの寄与はない
     */
//    AndoLab::Vector3d <double> rc2r1 = r1 - pmc.rc();
//    AndoLab::Vector3d <double> rc2r2 = r2 - pmc.rc();
//    if (
//        ( rc2r1.abs() > Coefficient_SDofPMC_accounted*pmc.sig_r() ) &&
//        ( rc2r2.abs() > Coefficient_SDofPMC_accounted*pmc.sig_r() ) &&
//        ( rc2r1%rc2r2 > 0.0 ) ){
//      delta += 0.0;
//    } else {
      /*
       * PMCの分布が含まれるのであれば、PMC用の区間で積分をしてゆく
//       */
//      Region pmc_region = search_pmc_region(r2, r1);

//      for(int i = 0; i < pmc_region.num; i++){

        /* 寄与を計算するベクトル */
//        AndoLab::Vector3d <double> length = pmc_region.r_to[i] - pmc_region.r_from[i];
//      AndoLab::Vector3d <double> length = r2 - r1;

        /* 分割点数。Interval毎 */
//        const int Num = int( length.abs() / PMC_Integral_Interval + 1 );

//        AndoLab::Vector3d <double> dr = length / double(Num); /* 分割距離 */
//        AndoLab::Vector3d <double> pre_rp = pmc_region.r_from[i];
//        AndoLab::Vector3d <double> pre_rp = r1;

//        for(int i_dr = 0; i_dr < Num; i_dr++){
//          AndoLab::Vector3d <double> r_p = pmc_region.r_from[i] + (i_dr+0.5)*dr;
//          AndoLab::Vector3d <double> r_p = r1 + (i_dr+0.5)*dr;

//          ofs_log << (r_p.abs() - Radius_of_Earth)*1e-3 << " "
//              << r_p.theta() << " " << pmc.dist( r_p ) << std::endl;

          /* r_p がPMC分布から近いときのみ計算する */
//          if ( (r_p - pmc.rc()).abs() < Coefficient_SDofPMC_accounted*pmc.sig_r() ){
//
//            if ( logging ){
//              ofs_log << pre_rp.x()*m2km << " "
//                  << pre_rp.y()*m2km << " "
//                  << pre_rp.z()*m2km << " "
//                  << (r_p-pre_rp).x()*m2km << " "
//                  << (r_p-pre_rp).y()*m2km << " "
//                  << (r_p-pre_rp).z()*m2km << " "
//                  << (pre_rp.abs() - Radius_of_Earth)*m2km << " "
//                  << (r_p.abs() - Radius_of_Earth)*m2km << " "
//                  << std::endl;
//            }

            //            delta += Mie_opt_depth(pre_rp, r_p, pmc);
//            delta += vector_gauss_quadrature(mie_beta, (void*)(&pmc), pre_rp, r_p);

//          }
//          pre_rp = r_p;
//        }

//        if ( (pre_rp - pmc.rc()).abs() < Coefficient_SDofPMC_accounted*pmc.sig_r() ){
//          delta += Mie_opt_depth(pre_rp, pmc_region.r_to[i], pmc); /* 最後の半区間の積分 */
//          delta += Mie_opt_depth(pre_rp, r2, pmc); /* 最後の半区間の積分 */
//          delta += vector_gauss_quadrature(mie_beta, (void*)(&pmc), pre_rp, r2); /* 最後の半区間の積分 */
    delta += vector_gauss_quadrature(mie_beta, (void*)(&pmc), r1, r2);
//        }

//      }
//    }
  }

  return delta;

}
