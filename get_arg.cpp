/*
 * get_arg.cpp
 *
 *  Created on: 2024/02/28
 *      Author: ando
 *
 * 引数から、解析する緯度・経度、日付・時刻を得る
 *
 * 引数の仕様
 * (実行ファイル) 緯度 経度 YYMMDD(or YYYYMMDD) hhmm LatLonIndex
 * argv[0]        [1]  [2]  [3]                [4]  [5]
 */


#include <iostream>
#include <string>

#include <Msis21.h>

#include "Date.h"
#include "pmc_simulation.h"

void get_arg(const int argc, char **argv, double &latitude, double &longitude,
    Date &date, AndoLab::Msis21 &msis, std::string &data_dir,
    std::string &process_id){

  if ( argc < 6 ){
    std::cerr << "Insufficient arguments." << std::endl;
    exit(1);
  }

  process_id = "";

  /* 緯度、経度 設定 */
  latitude = std::stod( std::string( argv[1] ) );
  longitude = std::stod( std::string( argv[2] ) );

  process_id += std::string( argv[1] ) + "_" + std::string( argv[2] ) + "_";

  std::string date_arg = argv[3];
  std::string YYMMDD;
  std::string YYYYMMDD;
  if ( date_arg.size() == 6 ){
    YYMMDD = date_arg;
    YYYYMMDD = "20" + YYMMDD;
  } else if ( date_arg.size() == 8 ){
    YYYYMMDD = date_arg;
    YYMMDD = YYYYMMDD.substr(2,6);
  } else {
    std::cerr << "Date must be YYMMDD or YYYYMMDD." << std::endl;
    exit(1);
  }
  std::string hhmm = argv[4];
  std::string YYYY = YYYYMMDD.substr(0,4);
  std::string YYYYMM = YYYYMMDD.substr(0,6);

  process_id += YYYYMMDD + "_" + hhmm + "_";

  /* データ保存ディレクトリ */
  data_dir = "data/" + YYYY + "/" + YYYYMM + "/" + YYYYMMDD;
  output_file_id = YYYYMMDD + "_" + hhmm + "_" + std::string( argv[5] );

  process_id += argv[5];

  int year = std::stoi( YYYYMMDD.substr(0,4) );
  int month = std::stoi( YYYYMMDD.substr(4,2) ); /* 月 */
  int day = std::stoi( YYYYMMDD.substr(6,2) );
  int hour = std::stoi( hhmm.substr(0,2) );
  int minute = std::stoi( hhmm.substr(2,2) );

  msis.date( year, month, day, hour*60 + minute );

  date.doy( msis.doy() );
  date.minute_day( hour*60 + minute );

  return;
}

//  std::string str_hour;
//  if ( argc > 1 ){
//    str_hour = argv[1];
//    Minute_of_Day = int( std::stod(str_hour) * 60 + 0.5 );
//  }
