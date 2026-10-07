/*
 * write_log.h
 *
 *  Created on: 2025/03/04
 *      Author: ando
 */

#ifndef WRITE_LOG_H_
#define WRITE_LOG_H_

#include <string>

extern std::ofstream ofs_log; /* ログ出力 */
extern LOGGING logging_method; /* ログ出力の方法 */
extern std::string logging_str;


void write_log_open(std::string data_directory);
void write_log_close(std::string data_directory);


#endif /* WRITE_LOG_H_ */
