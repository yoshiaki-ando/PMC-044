/*
 * region.h
 *
 *  Created on: 2019/12/27
 *      Author: ando
 */

#ifndef REGION_H_
#define REGION_H_

#include <Vector3d.h>

class Region{
public:
  Region(void): num { 0 } { };
  int num { 0 };
  AndoLab::Vector3d <double> r_from[2];
  AndoLab::Vector3d <double> r_to[2];
  void push(AndoLab::Vector3d <double> vr_from, AndoLab::Vector3d <double> vr_to){
    r_from[num] = vr_from;
    r_to[num] = vr_to;
    num++;
  }
};



#endif /* REGION_H_ */
