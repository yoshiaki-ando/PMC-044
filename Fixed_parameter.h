/*
 * Fixed_parameter.h
 *
 *  Created on: 2024/05/27
 *      Author: ando
 */
#include <cmath>


#ifndef FIXED_PARAMETER_H_
#define FIXED_PARAMETER_H_

constexpr double m2km { 1.0e-3 };

constexpr double Deg2rad { M_PI / 180.0 };
constexpr double Rad2deg { 180.0 / M_PI };

constexpr int IDX_LATITUDE { 0 };
constexpr int IDX_LONGITUDE { 1 };

constexpr double Radius_of_Earth { 6370.e3 };
constexpr double Altitude_of_GEO { 35790.e3 };
constexpr double Rgeo { Altitude_of_GEO + Radius_of_Earth };

constexpr double Longitude_of_himawari { 140.7 };

/* 使用する波長 */
constexpr int Num_Lambda { 3 };
constexpr double Lambda[Num_Lambda] { 470e-9, 510e-9, 640e-9 }; /* 470nm, 510nm, 640nm */

constexpr double Lower_PMC { 78.0e3 };
constexpr double Upper_PMC { 88.0e3 };
constexpr double DAlt { 1.0e3 };
constexpr int idx_LowerPMC { int(Lower_PMC/DAlt) };
constexpr int idx_UpperPMC { int(Upper_PMC/DAlt) };
constexpr int Num_PMC_Layer { idx_UpperPMC - idx_LowerPMC + 1 };

constexpr double Altitude_of_Atmosphere { 100.e3 };



constexpr double R_lower_pmc { Radius_of_Earth + Lower_PMC };
constexpr double R_upper_pmc { Radius_of_Earth + Upper_PMC };

/* 光学的深さを求める際の、領域を分別する記号 */
constexpr bool PMC_REGION { true };
constexpr bool AIR_REGION { false };

/* LOS上の寄与を計算する区間がどこにあるか */
constexpr int HIGHER_THAN_PML_LAYER { 0 };
constexpr int INSIDE_PML_LAYER { 1 };
constexpr int LOWER_THAN_PML_LAYER { 2 };

/* LOSのδを求めるのか、太陽方向のδを求めるのか */
//constexpr

#endif /* FIXED_PARAMETER_H_ */
