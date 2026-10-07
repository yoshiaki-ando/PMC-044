/*
 * select_latlon.cpp
 *
 *  Created on: 2024/06/11
 *      Author: ando
 */
#include <iostream>
#include <cmath>

#include "Fixed_parameter.h"
#include "pmc_observation.h"

int select_latlon(double **latlon, const double Latitude, const double Longitude){
  int idx = -1; /* 探した latlon のインデックス */
  double min_err = 1.0e8; /* 指定した経度と、データの経度との最小の誤差のものを返す */

  for(int i = 0; i < Num_observed_data_latitude; i++){

    if ( std::abs( latlon[IDX_LATITUDE][i] - Latitude ) < 1e-3 ){

      double err = std::abs( latlon[IDX_LONGITUDE][i] - Longitude );

      if ( err < min_err ){
        min_err = err;
        idx = i;
      }

    }
  }

  if ( idx == -1 ){
    std::cout << "Lat, Lon (" << Latitude << ", " << Longitude << " ) is not found." << std::endl;
    exit(2);
  }

  return idx;
}

