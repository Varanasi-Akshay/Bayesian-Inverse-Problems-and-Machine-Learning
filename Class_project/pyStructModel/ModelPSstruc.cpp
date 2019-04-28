//
//  ModelPSstruc.cpp
//  fullstructuralmodel
//
//  Created by Will Zhang on 6/12/18.
//  Copyright © 2018 Will Zhang. All rights reserved.
//

#include "ModelPSstruc.hpp"
#include <math.h>



void model_ps_struc::calc_new_old_fraction(double rate, double time)
{
    rate_PS = rate;
    time_elapsed = time;
    
    fraction_old = exp( - rate * time);
    fraction_new = 1.0 - fraction_old;
}

void model_ps_struc::update_reference(double arg[4])
{
    F11_pushforward = arg[0];
    F12_pushforward = arg[1];
    F21_pushforward = arg[2];
    F22_pushforward = arg[3];
}

void model_ps_struc::update_history(double arg[4])
{
    for (int i = 0; i < 4; i++) {
        strainhistory[i] = arg[i];
    }
}

void model_ps_struc::updata_parameters(double arg[13])
{
    kM = arg[0];
    ratio = arg[1];
    alpham1 = arg[2];
    betam1 = arg[3];
    kC = arg[4];
    ODFmean = arg[5];
    ODFstdev = arg[6];
    Dxmean = arg[7];
    Dxstdev = arg[8];
    Dxlb = arg[9];
    Dxub = arg[10];
    Dx_aniso = arg[11];
    kI = arg[12];
}


model_ps_struc::model_ps_struc(double arg[]):
    matrixmodel(&arg[0]), collagenmodel(&arg[4]),
    kM{arg[0]}, ratio{arg[1]}, alpham1{arg[2]}, betam1{arg[3]},
    kC{arg[4]}, ODFmean{arg[5]}, ODFstdev{arg[6]},
    Dxmean{arg[7]}, Dxstdev{arg[8]}, Dxlb{arg[9]}, Dxub{arg[10]}, Dx_aniso{arg[11]},
    kI{arg[12]},
    F11_pushforward{arg[13]}, F12_pushforward{arg[14]},F21_pushforward{arg[15]},F22_pushforward{arg[16]},
    strainhistory{arg[19], arg[20], arg[21], arg[22]}
    {calc_new_old_fraction(arg[17], arg[18]);};




int model_ps_struc::stress(double vF[4], double res[4], bool crosslinked, bool withoutmatrix)
{
    double identityM[4] = {1.0, 0.0, 0.0, 1.0};
    
    double collagenstress[4];
    double oldmatrixstress[4];
    double newmatrixstress[4];
    
    
    matrixmodel.stress(identityM, vF, oldmatrixstress);
    matrixmodel.stress(strainhistory, vF, newmatrixstress);
    
    if (crosslinked) {
        // Using the full exogenously crosslinked structural model
        // including quadruple integral of interaction term
        collagenmodel.stress_total(vF, collagenstress);
    } else {
        // Using only the standard structural model
        // without quadruple integral of interaction term
        collagenmodel.stress(vF, collagenstress);
    }
    
    for (int i = 0; i < 4; i++) {
        if (withoutmatrix){
            res[i] = collagenstress[i];
        } else {
            res[i] = collagenstress[i] + fraction_old * oldmatrixstress[i] + fraction_new * newmatrixstress[i];
        }
        
    }
    
    return 0;
}
















