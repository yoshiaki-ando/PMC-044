/*
 * across.h
 *
 *  Created on: 2019/12/02
 *      Author: ando
 */

#ifndef ACROSS_PTS_H_
#define ACROSS_PTS_H_

#include <Vector3d.h>

class Across_pts{
private:
public:
  Across_pts(AndoLab::Vector3d <double> r2, AndoLab::Vector3d <double> r1,
      double Altitude);
  AndoLab::Vector3d <double> r[2];
  int num;

  void copy( /* r をコピーする */
      AndoLab::Vector3d <double> *destination /* コピー先 */
  );
};



#endif /* ACROSS_PTS_H_ */
