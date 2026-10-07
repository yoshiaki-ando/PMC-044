/*
 * write_log.cpp
 *
 *  Created on: 2025/03/04
 *      Author: ando
 */
#include <iostream>
#include <fstream>
#include <string>

#include "pmc_simulation.h"
#include "write_log.h"

namespace{
std::string log_filename(){
  if ( output_file_id.empty() ){
    return "log.txt";
  }
  return "log_" + output_file_id + ".txt";
}
}

void write_log_open(std::string data_dir){
  std::string str_output = "# Iteration Error CenterAlt[km] Sig_z[km] Sig_r[km] R0[nm] RadiusWidth[nm] Nc[m^-2] N0[x10^6 m^-3] Shift[km]\n";

  if ( logging_method == LOGGING::COLLECTED ){
    logging_str += str_output;
  } else if ( logging_method == LOGGING::SUCCESSIVE ){
    ofs_log.open( data_dir + "/" + log_filename(), std::ios::out);
    ofs_log << "# Iteration Error CenterAlt[km] Sig_z[km] Sig_r[km] R0[nm] RadiusWidth[nm] Nc[m^-2] N0[x10^6 m^-3] Shift[km]\n";
  }

}


void write_log_close(std::string data_dir){

  if ( logging_method == LOGGING::COLLECTED ){
    ofs_log.open( data_dir + "/" + log_filename(), std::ios::out);
    ofs_log << logging_str;
  }

  if ( logging_method != LOGGING::NO ){
    ofs_log.close();
  }
}
