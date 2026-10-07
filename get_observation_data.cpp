/*
 * get_observation_data.cpp
 *
 *  Created on: 2024/06/07
 *      Author: ando
 */
/*
 * 観測データを読み込む
 */
#include <iostream>
#include <fstream>
#include <string>

#include <Command.h>

#include "Fixed_parameter.h"
#include "pmc_observation.h"

void get_observation_data(
    double ***intensity,
    double **sigma,
    double **latlon,
    char **argv
    ){

  std::string line;
  std::string date_arg( argv[3] );
  std::string YYYYMMDD;
  if ( date_arg.size() == 6 ){
    YYYYMMDD = "20" + date_arg;
  } else if ( date_arg.size() == 8 ){
    YYYYMMDD = date_arg;
  } else {
    std::cerr << "Date must be YYMMDD or YYYYMMDD." << std::endl;
    exit(1);
  }
  std::string hhmm( argv[4] );

  constexpr int Offset_for_not_altitude_data { 2 }; /* 最初の2行は緯度・経度情報 */

  for(int j_lambda = 0; j_lambda < Num_Lambda; j_lambda++){
//    std::ifstream ifs( ("data/h08_b0" + std::to_string(j_lambda+1) + "_s01s02_20160709_210000.txt").c_str() );
    std::ifstream ifs( (obs_data_dir + "/h08_b0" + std::to_string(j_lambda+1) +
        "_s01s02_" + YYYYMMDD + "_" + hhmm + "00.txt").c_str() );
    if ( !ifs ){
      std::cout << "File: \""
          << obs_data_dir + "/h08_b0" + std::to_string(j_lambda+1) +
          "_s01s02_" + YYYYMMDD + "_" + hhmm + "00.txt\" does not exist." << std::endl;
      exit(-1);
    }

    std::getline(ifs, line); /* 1行目はただ height */
    std::getline(ifs, line); /* 2行目は lat lon 0 〜 100 */

    for(int i_alpha = 0; i_alpha < Num_observed_data_latitude; i_alpha++){
      std::getline(ifs, line); /* 3行目以降は lat lon 0 〜 100の輝度 */
      AndoLab::Command cmd(line, ' ');

      /* 緯度・経度情報を記憶 */
      latlon[IDX_LATITUDE][i_alpha] = cmd.get_d(0);
      latlon[IDX_LONGITUDE][i_alpha] = cmd.get_d(1);

      for(int k = 0; k < Num_observed_data_altitude; k++){
        intensity[j_lambda][i_alpha][k] = cmd.get_d(Offset_for_not_altitude_data + k);
      }
    }
    ifs.close();
  }

  /* 標準偏差 */

  int j_lambda = 0;
  std::string line_sigma;
  std::ifstream ifs_sigma( (obs_data_dir + "/sigma_b0" + std::to_string(j_lambda+1) +
      "_s01s02_" + YYYYMMDD + "_" + hhmm + "00.txt").c_str() );

  if ( !ifs_sigma ){
    std::cout << "File: \""
        << obs_data_dir + "/sigma_b0" + std::to_string(j_lambda+1) +
        "_s01s02_" + YYYYMMDD + "_" + hhmm + "00.txt\" does not exist."
        << std::endl;
    exit(-1);
  }

  std::getline(ifs_sigma, line_sigma); /* 1行目はただ height */
  std::getline(ifs_sigma, line_sigma); /* 2行目は lat lon 0 〜 100 */
  for(int i_alpha = 0; i_alpha < Num_observed_data_latitude; i_alpha++){
    std::getline(ifs_sigma, line_sigma); /* 3行目以降は lat lon 0 〜 100の標準偏差 */
    AndoLab::Command cmd_sigma(line_sigma, ' ');
    for(int k = 0; k < Num_observed_data_altitude; k++){
      sigma[i_alpha][k] = cmd_sigma.get_d(Offset_for_not_altitude_data + k);
    }
  }
  ifs_sigma.close();


}

