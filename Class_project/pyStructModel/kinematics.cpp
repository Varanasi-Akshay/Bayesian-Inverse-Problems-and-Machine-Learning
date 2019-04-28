//
//  kinematics.cpp
//  effectivemodel
//
//  Created by Will Zhang on 12/12/17.
//  Copyright © 2017 Will Zhang. All rights reserved.
//

#include <math.h>
#include "kinematics.hpp"


//  The dot product between a matrix and a vector
void dot2D(double vF[4], double vM[2], double res[2]) {
    
    for (int i = 0; i < 2; i++) {
        res[i] = vF[2*i] * vM[0] + vF[2*i+1] * vM[1];
    }
    
    return;
}




//  multiplying 2x2 matrix together
void mult2D(double vT1[4], double vT2[4], double res[4]) {
    
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            res[2*i+j] = 0.0;
        }
    }
    
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                res[2*i+j] += vT1[2*i+k]*vT2[2*k+j];
            }
        }
    }
    
    return;
}




void mult2Dv(double vT1[4], double vM[2], double res[2]) {
    
    for (int i = 0; i < 2; i++) {
        res[i] = 0.0;
    }
    
    for (int i = 0; i < 2; i++) {
        for (int k = 0; k < 2; k++) {
            res[i] +=vT1[2*i+k]*vM[k];
        }
    }
    
    return;
}





void cross2D(double vM[2], double res[2]){
    res[0] = -vM[1];
    res[1] = vM[0];
    return;
}





void trans2D(double vT1[4], double res[4]) {
    res[0] = vT1[0];
    res[1] = vT1[2];
    res[2] = vT1[1];
    res[3] = vT1[3];
    return;
}



void tensorproduct2D(double v1[2], double v2[2], double res[4]) {
    
    for (int i = 0; i < 2 ; i++) {
        for (int j = 0; j < 2; j++) {
            res[2*i+j] = v1[i] * v2[j];
        }
    }
    
    return;
}






void calc_m(double vF[4], double vM[2], double res[2]) {
    
    for (int i = 0; i < 2; i++) {
        res[i] = 0.0;
    }
    
    for (int i = 0; i < 2; i++) {
        for (int k = 0; k < 2; k++) {
            res[i] += vF[2*i+k]*vM[k];
        }
    }

    double norm = 0.0;
    
    for (int i = 0; i < 2; i++) {
        norm += res[i]*res[i];
    }
    
    norm = sqrt(norm);
    
    for (int i = 0; i < 2; i++) {
        res[i] = res[i] / norm;
    }
    
    return;
}





void calc_s(double vF[4], double vM[2], double res[2]){
    
    double vm[2];
    calc_m(vF, vM, vm);
    cross2D(vm, res);
    
    return;
}





double calc_lambda_M(double vF[4], double vM[2]) {
    double res = 0.0;
    
    double vm[2];
    
    calc_m(vF, vM, vm);
    
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            res += vm[i] * vF[2*i+j] * vM[j];
        }
    }
    
    return res;
}




double calc_lambda_S(double vF[4], double vM[2]) {
    double res = 0.0;
    
    double vs[2], vS[2];
    
    cross2D(vM, vS);
    calc_s(vF, vM, vs);
    
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            res += vs[i] * vF[2*i+j] * vS[j];
        }
    }
    
    return res;
}




double calc_phi(double vF[4], double vM[2]){
    double res = 0.0;
    
    double vm[2], vS[2];
    
    double lambda_M;
    
    lambda_M = calc_lambda_M(vF, vM);
    
    cross2D(vM, vS);
    calc_m(vF, vM, vm);
    
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            res += vm[i] * vF[2*i+j] * vS[j];
        }
    }
    
    return res/lambda_M;
}




void calc_y(double vF[4], double vM[2], double res[3]) {
    double vm[2], vs[2], vS[2];
    
    cross2D(vM, vS);
    calc_m(vF, vM, vm);
    cross2D(vm, vs);
    
    for (int i = 0; i < 3; i++) {
        res[i] = 0.0;
    }
    
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            res[0] += vm[i] * vF[2*i+j] * vM[j];
            res[1] += vs[i] * vF[2*i+j] * vS[j];
            res[2] += vm[i] * vF[2*i+j] * vS[j];
        }
    }
    
    res[2] = res[2]/res[0];
    res[1] = log(res[1]);
    res[0] = log(res[0]);
    
    return;
}



void calc_E(double vF[4], double vM[2], double res[3]) {
    
    double vS[2];
    cross2D(vM, vS);
    
    for (int i = 0; i < 3; i++) {
        res[i] = 0.0;
    }
    
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                res[0] += vM[i] * vF[2*j+i] * vF[2*j+k] * vM[k];
                res[1] += vS[i] * vF[2*j+i] * vF[2*j+k] * vS[k];
                res[2] += vM[i] * vF[2*j+i] * vF[2*j+k] * vS[k];
            }
        }
    }
    
    res[0] = 0.5*(res[0] - 1.0);
    res[1] = 0.5*(res[1] - 1.0);
    res[2] = 0.5*(res[2]);
    
    
    return;
}



void calc_C_from_F(double vT[4], double res[4])
{
    double vTinv[4];
    
    trans2D(vT, vTinv);
    
    mult2D(vTinv, vT, res);
}


void calc_B_from_F(double vT[4], double res[4])
{
    double vTinv[4];
    
    trans2D(vT, vTinv);
    
    mult2D(vT, vTinv, res);
}


void basis_y_T(double vF[4], double vM[2], double res[3][4]) {
    double vm[2], vs[2], vS[2];
    calc_m(vF, vM, vm);
    calc_s(vF, vM, vs);
    cross2D(vM, vS);
    
    double lambda_M = 0.0, lambda_S = 0.0;
    
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            lambda_M += vm[i] * vF[2*i+j] * vM[j];
            lambda_S += vs[i] * vF[2*i+j] * vS[j];
        }
    }
    
    for (int i = 0; i < 2 ; i++) {
        for (int j = 0; j < 2; j++) {
            res[0][2*i+j] = vm[i] * vm[j];
            res[1][2*i+j] = vs[i] * vs[j];
            res[2][2*i+j] = (lambda_M/lambda_S) * ( vm[i] * vs[j] + vs[i] * vm[j] );
        }
    }
    
    return;
}








void basis_y_S(double vF[4], double vM[2], double res[3][4]) {
    
    double vm[2], vs[2], vS[2];
    cross2D(vM, vS);
    
    double lambda_M = 0.0, lambda_S = 0.0, phi = 0.0;
    
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            lambda_M += vm[i] * vF[2*i+j] * vM[j];
            lambda_S += vs[i] * vF[2*i+j] * vS[j];
            phi += vm[i] * vF[2*i+j] * vS[j];
        }
    }
    
    double lambda_M2;
    
    lambda_M2 = lambda_M*lambda_M;
    
    for (int i = 0; i < 2 ; i++) {
        for (int j = 0; j < 2; j++) {
            res[0][2*i+j] = vM[i] * vM[j] / lambda_M2;
            res[1][2*i+j] = (vS[i] * vS[j] - phi * (vM[i] * vS[j] + vS[i] * vM[j]) + phi*phi* vM[i] * vM[j]) / (lambda_S * lambda_S);
            res[2][2*i+j] = ( vM[i] * vS[j] + vS[i] * vM[j] - 2.0 * phi * vM[i] * vM[j] ) / lambda_M2;
        }
    }
    
    
    return;
}





void basis_E_T(double vF[4], double vM[2], double res[3][4]) {
    
    double vm[2] = {0.0, 0.0}, vs[2] = {0.0, 0.0}, vS[2];
    cross2D(vM, vS);
    
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            vm[i] += vF[2*i+j]*vM[j];
            vs[i] += vF[2*i+j]*vS[j];
        }
    }
    
    for (int i = 0; i < 2 ; i++) {
        for (int j = 0; j < 2; j++) {
            res[0][2*i+j] = vm[i] * vm[j];
            res[1][2*i+j] = vs[i] * vs[j];
            res[2][2*i+j] = 0.5 * (vm[i]*vs[j]  +  vs[i]*vm[j]);
        }
    }
    
    
    return;
}






void basis_E_S(double vF[4], double vM[2], double res[3][4]) {
    
    double vS[2];
    cross2D(vM, vS);
    
    for (int i = 0; i < 2 ; i++) {
        for (int j = 0; j < 2; j++) {
            res[0][2*i+j] = vM[i] * vM[j];
            res[1][2*i+j] = vS[i] * vS[j];
            res[2][2*i+j] = 0.5 * (vM[i]*vS[j]  +  vS[i]*vM[j]);
        }
    }
    
    return;
}






