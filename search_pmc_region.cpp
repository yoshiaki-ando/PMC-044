#include <iostream>

#include <Vector3d.h>

#include "pmc_simulation.h"
#include "Region.h"
#include "Across_pts.h"

Region search_pmc_region(
    AndoLab::Vector3d <double> &r2,
    AndoLab::Vector3d <double> &r1
    ){

  Region v;

  Across_pts across_lower_pmc(r2, r1, Lower_PMC);
  Across_pts across_upper_pmc(r2, r1, Upper_PMC);
  const int num = across_lower_pmc.num + across_upper_pmc.num;

  if ( num == 0 ){
    v.num = 0;
    return v;
  }

  /* PMC領域の外側 */
  if ( (r2.r() < R_lower_pmc) || (r2.r() > R_upper_pmc) ){

    /* ケース4: PMC領域の外側で、r1はPMC領域内
     PMC領域へ内側から入る */
    if ( (num == 1) && (across_lower_pmc.num == 1) ){
//      std::cout << "case a: ";
      v.push(across_lower_pmc.r[0], r1);

      /* PMC領域へ外側から入る */
    } else if ( num == 1 ){
//      std::cout << "case b: ";
      v.push(across_upper_pmc.r[0], r1);

      /* ケース2: PMC領域の外側で、r2->r1でPMC領域を1つ抜ける */
      /* 内←→外と抜ける場合 */
    } else if ( (num == 2) && (across_lower_pmc.num == 1) ){
//      std::cout << "case c: ";
      v.push(across_lower_pmc.r[0], across_upper_pmc.r[0]);
//      delta_mie = Mie_opt_depth(across_lower_pmc.r[0], across_upper_pmc.r[0], pmc);

      /* 外→外と抜ける場合 */
    } else if ( num == 2 ){
//      std::cout << "case d: ";
      v.push(across_upper_pmc.r[0], across_upper_pmc.r[1]);
//      delta_mie = Mie_opt_depth(across_upper_pmc.r[0], across_upper_pmc.r[1], pmc);

      /* ケース5: PMC領域の外側で、PMC領域を一つ抜け、r1はPMC領域内 */
      /* PMCへは内側から入る */
    } else if ( (num == 3) && (across_lower_pmc.num == 2) ){
//      std::cout << "case f\n";
      v.push(across_lower_pmc.r[0], across_upper_pmc.r[0]);
      v.push(across_lower_pmc.r[1], r1);
//      delta_mie
//      = Mie_opt_depth(across_lower_pmc.r[0], across_upper_pmc.r[0], pmc) +
//      Mie_opt_depth(across_lower_pmc.r[1], r1, pmc);

      /* PMCへは外側から入る */
    } else if ( num == 3 ) {
      std::cout << "case g (impossible)\n";
      v.push(across_lower_pmc.r[0], across_upper_pmc.r[0]);
      v.push(across_upper_pmc.r[1], r1);
//      delta_mie
//      = Mie_opt_depth(across_lower_pmc.r[0], across_upper_pmc.r[0], pmc) +
//      Mie_opt_depth(across_upper_pmc.r[1], r1, pmc);

      /* ケース3: PMC領域の外側で、r2->r1でPMC領域が2つ抜ける */
    } else if ( num == 4 ) {
//      std::cout << "case h\n";
//      delta_mie
//      = Mie_opt_depth(across_lower_pmc.r[0], across_upper_pmc.r[0], pmc) +
//      Mie_opt_depth(across_lower_pmc.r[1], across_upper_pmc.r[1], pmc);
      v.push(across_lower_pmc.r[0], across_upper_pmc.r[0]);
      v.push(across_lower_pmc.r[1], across_upper_pmc.r[1]);
    } else {
      throw( std::string("(search_pmc_region) undetermined case.") );
    }


    /* PMC領域の内側 */
  } else {
    /* ケース6: PMC領域の内側で、r2->r1でPMC領域を抜ける */
    /* 内側へ抜ける */
    if ( (num == 1) && (across_lower_pmc.num == 1) ){
//      std::cout << "case i: ";
//      delta_mie = Mie_opt_depth(r2, across_lower_pmc.r[0], pmc);
      v.push(r2, across_lower_pmc.r[0]);

      /* 外側へ抜ける */
    } else if ( num == 1) {
//      std::cout << "case j: ";
//      delta_mie = Mie_opt_depth(r2, across_upper_pmc.r[0], pmc);
      v.push(r2, across_upper_pmc.r[0]);

      /* ケース7: PMC領域の内側で、r2->r1でPMC領域を抜けて、r1がPMC領域
       * upper=1 and lower=1 */
    } else if ( num == 2 ){
      /* PMC領域の内側へ抜けて、内側から入るケースしかありえない */
//      std::cout << "case k : Lower = " << across_lower_pmc.num
//          << " & Upper = " << across_upper_pmc.num
//          << " , (Lower must be 2)\n";
//      delta_mie
//      = Mie_opt_depth(r2, across_lower_pmc.r[0], pmc) +
//      Mie_opt_depth(across_lower_pmc.r[1], r1, pmc);
      v.push(r2, across_lower_pmc.r[0]);
      v.push(across_lower_pmc.r[1], r1);

    } else if ( num == 3 ){


      /* ケース8: PMC領域の内側で、r2->r1でPMC領域を抜けて、PMC領域をもう一つ
       * 抜ける
       * upper=1 and lower=2 */
//      std::cout << "case l\n";
//      delta_mie
//      = Mie_opt_depth(r2, across_lower_pmc.r[0], pmc)
//      + Mie_opt_depth(across_lower_pmc.r[1], across_upper_pmc.r[0], pmc);
      v.push(r2, across_lower_pmc.r[0]);
      v.push(across_lower_pmc.r[1], across_upper_pmc.r[0]);
    } else {
      throw( std::string("(search_pmc_region) undetermined case.") );
    }
  }

  /* ケース1: PMC領域の外側で、r2->r1でPMC領域がない
   * upper=0, lower=0 */

  return v;
}
