//
//  ModelCollagenFiber.cpp
//  fullstructuralmodel
//
//  Created by Will Zhang on 5/10/18.
//  Copyright © 2018 Will Zhang. All rights reserved.
//


#include "math.h"
#include "ModelCollagenFiber.hpp"



//Fiber Model that's linear in PF
double model_col_fiber_PF(double par[], int parn, double arg[2], int argn)
{
    // arg[0] = lambda_f, stretch
    // arg[1] = lambda_s, slack
    // Number of arguments is argn = 2
    
    
    double lf = 1.0 / arg[0];
    double ls = 1.0 / arg[1];
    
    return ls * (ls - lf);
}





//Fiber Model that's linear in SE
double model_col_fiber_SE (double par[], int parn, double arg[2], int argn)
{
    // arg[0] = lambda_f, stretch
    // arg[1] = lambda_s, slack
    // Number of arguments is argn = 2
    
    double Ef = arg[0]*arg[0];
    double Es = arg[1]*arg[1];
    double denom = 1.0/Es;
    
    return 0.5*(Ef - Es) * denom * denom;
};









