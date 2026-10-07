/*
 * Date.h
 *
 *  Created on: 2019/07/09
 *      Author: ando
 */

#ifndef DATE_H_
#define DATE_H_

class Date{
private:
  int pDoy, /* doy */
  pMinute;  /* 0 - (24*60-1) */
public:
  Date(void){ };
  Date(int doy, int minute_in_a_day):
    pDoy(doy), pMinute(minute_in_a_day){ };

  void doy(int v_doy){ pDoy = v_doy; }
  int doy(void){ return pDoy; }
  void minute_day(int v_minute){ pMinute = v_minute; }
  int minute_day(void){ return pMinute; }

  /* 夏至からのDOY。夏至は常に6/21と仮定 */
  int doy_from_solstice(void){
    int doy = pDoy - 172; /* 172...6/21 */
    if ( doy < 0 ){
      doy += 365;
    }
    return doy;
  }
};



#endif /* DATE_H_ */
