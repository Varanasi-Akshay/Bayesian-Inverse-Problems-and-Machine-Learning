//
//  ModelCollagenTissueExp.cpp
//  fullstructuralmodel
//
//  Created by Will Zhang on 5/11/18.
//  Copyright © 2018 Will Zhang. All rights reserved.
//



#include "ModelCollagenTissueExp.hpp"


/*****          Exponential type ensemble models           *****/

double col_ens_PF_exp::stress(double lf)
{
    if (lf>1.0) {
        return (exp(parB * 0.5*(lf*lf - lmax)) - scaling);
    } else {
        return 0.0;
    }
}

double col_ens_PF_exp::strainenergy(double lf)
{
    double lfm1 = 0.5*(lf*lf - 1.0);
    if (lf > 1.0) {
        return scaling*( (exp(parB * lfm1) - 1.0)/parB - lfm1);
    } else {
        return 0.0;
    }
}
