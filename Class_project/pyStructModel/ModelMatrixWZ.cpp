//
//  ModelMatrixWZ.cpp
//  fullstructuralmodel
//
//  Created by Will Zhang on 5/14/18.
//  Copyright © 2018 Will Zhang. All rights reserved.
//

#include "math.h"
#include "kinematics.hpp"

#include "ModelMatrixWZ.hpp"

//  The default constructor
model_matrix_yeoh_wz::model_matrix_yeoh_wz():
kM{1.0}, ratio{1.0}, alpham1{1.0}, betam1{1.0}
{
}




//  Parameterized constructor
model_matrix_yeoh_wz::model_matrix_yeoh_wz(double arg[4]):
//  Below are the parameters of the model
kM{arg[0]}, //  This is the modulus
ratio{arg[1]},  //  The is the ratio between the linear and the non-linear part
alpham1{arg[2] - 1.0},    //  This is the power of the linear part
betam1{arg[3] - 1.0}      //  This is the power of the non-linear part
{
}


void model_matrix_yeoh_wz::calc_tensorinverse(double vT[4], double vTinv[4], double &det)
{
    det = vT[0]*vT[3] - vT[1]*vT[2];
    
    vTinv[0] = vT[3]/det;
    vTinv[1] = -vT[1]/det;
    vTinv[2] = -vT[2]/det;
    vTinv[3] = vT[0]/det;
}


double model_matrix_yeoh_wz::calc_invariant_1(double vReferencestate[4], double vC[4])
{
    double res;
    
    //  compute the trace of (A^-T C A^-1) = A_j1 * C_jk * A_k1 + A_j2 * C_jk * A_k2
    //  where A is the strain history, or the pushforward tensor vFinv_pushforward
    res = 0;
    for (int j = 1; j < 4; j++) {
        for (int k = 1; k < 4; k++) {
            res += (vReferencestate[2*j + 0] * vC[2*j + k] * vReferencestate[2*k + 0]
                    + vReferencestate[2*j + 1] * vC[2*j + k] * vReferencestate[2*k + 1]);
        }
    }
    return 0.0;
}


//  Compute the stresses
int model_matrix_yeoh_wz::stress(double vReferencestateA[4], double vF[4], double res[4])
{
    
    double vAinv[4];
    double detA = 1.0;
    double vBinv[4];
    double detB = 1.0;
    double vC[4];
    double vCinv[4];
    double detC = 1.0;
    double C33;
    double I1m3;
    double dPsi_dI;
    
    double dummieConstant;
    
    this->calc_tensorinverse(vReferencestateA, vAinv, detA);  //  Calculate the inverse of the strain history
    
    calc_C_from_F(vAinv, vBinv);    //  Note we are dealing with inverses, thus Binv to Finv^T Finv
    
    detB = detA * detA;    //  Get the determinant of Binv, which is 1/Binv_33
    
    calc_C_from_F(vF, vC);
    
    this->calc_tensorinverse(vC, vCinv, detC);
    
    C33 = 1.0 / detC;
    
    I1m3 = calc_invariant_1(vC, vC) + C33 - 3.0;
    
    dPsi_dI = kM * (pow(I1m3, alpham1) + ratio * pow(I1m3, betam1));
    
    dummieConstant = C33 * detB;   //  Note here that detBinv = 1.0/Binv_33
    
    res[0] = dPsi_dI * (vBinv[0] - dummieConstant * vCinv[0]);
    res[1] = dPsi_dI * (vBinv[1] - dummieConstant * vCinv[1]);
    res[2] = dPsi_dI * (vBinv[2] - dummieConstant * vCinv[2]);
    res[3] = dPsi_dI * (vBinv[3] - dummieConstant * vCinv[3]);
    
    return 0;
}













