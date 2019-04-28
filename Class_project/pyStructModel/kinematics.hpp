//
//  kinematics.hpp
//  effectivemodel
//
//  Created by Will Zhang on 12/12/17.
//  Copyright © 2017 Will Zhang. All rights reserved.
//

#ifndef kinematics_hpp
#define kinematics_hpp

void dot2D(double vF[4], double vM[2], double res[2]);

void mult2D(double vT1[4], double vT2[4], double res[4]);

void cross2D(double vM[2], double res[2]);

void trans2D(double vT1[4], double res[4]);

void tensorproduct2D(double v1[2], double v2[2], double res[4]);

void calc_m(double vF[4], double vM[2], double res[2]);

void calc_s(double vF[4], double vM[2], double res[2]);

double calc_lambda_M(double vF[4], double vM[2]);

double calc_lambda_S(double vF[4], double vM[2]);

double calc_phi(double vF[4], double vM[2]);

void calc_y(double vF[4], double vM[2], double res[3]);

void calc_E(double vF[4], double vM[2], double res[3]);

void calc_C_from_F(double vT[4], double res[4]);

void calc_C_from_F(double vT[4], double res[4], double det);

void calc_B_from_F(double vT[4], double res[4]);

void basis_y_T(double vF[4], double vM[2], double res[3][4]);

void basis_y_S(double vF[4], double vM[2], double res[3][4]);

void basis_E_T(double vF[4], double vM[2], double res[3][4]);

void basis_E_S(double vF[4], double vM[2], double res[3][4]);


#endif /* kinematics_hpp */
