/*
 * this_position_is_near_pmc.cpp
 *
 *  Created on: 2024/06/18
 *      Author: ando
 *
 * r1 → r2 の経路中、pmc の分布と近くなる点が存在するか
 */
#include "pmc_simulation.h"

/* PMCの光学的深さを考慮するために、PMC中心からPMC分布の何SDまで計算に入れるか
 * これを超えた経路については計算をしない */
constexpr double Coefficient_SDofPMC_accounted { 3.0 };

double get_nearest_point(
    AndoLab::Vector3d <double> rc,
    AndoLab::Vector3d <double> &r1,
    AndoLab::Vector3d <double> &r2){

  AndoLab::Vector3d <double> Rc = rc - r1;
  AndoLab::Vector3d <double> R2 = r2 - r1;

  const double sR2 = R2.abs();

  return Rc%R2 / sR2 / sR2;
}

bool This_position_is_near_PMC(PMC &pmc, AndoLab::Vector3d <double> &r1, AndoLab::Vector3d <double> &r2 ){

  bool is_near { false };
  const double t = get_nearest_point( pmc.rc(), r1, r2 );

  AndoLab::Vector3d <double> Rc = pmc.rc() - r1;
  double distance;
  if ( t <= 0.0 ){
    distance = Rc.abs();

  } else if ( t < 1.0 ){
    double sRc = Rc.abs();
    double Rc_nR2 = Rc%( (r2 - r1).n() );
    distance = std::sqrt( sRc*sRc - Rc_nR2*Rc_nR2 );

  } else {
    distance = ( Rc - (r2 - r1) ).abs();
  }

//  std::cout << "t = " << t << " : dist = " << distance << " ? "
//      << Coefficient_SDofPMC_accounted * pmc.sig_r() << std::endl;

  if ( distance < Coefficient_SDofPMC_accounted * pmc.sig_r() ){
    is_near = true;
  }

  return is_near;
}


