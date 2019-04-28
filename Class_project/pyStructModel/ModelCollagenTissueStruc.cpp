//
//  ModelCollagenTissueStruc.cpp
//  fullstructuralmodel
//
//  Created by Will Zhang on 5/11/18.
//  Copyright © 2018 Will Zhang. All rights reserved.
//


#include "kinematics.hpp"

#include "ModelCollagenTissueStruc.hpp"



int model_col_tissue_struc_PF::stress(double vF[4], double res[4])
{
    
    double vC[4];
    
    double lambda_ens[31];
    double Sf[31]; //Ensemble stress at each theta
    
    calc_C_from_F(vF, vC);
    
    for (int i = 0; i < 31; i++) {
        lambda_ens[i] = sqrt(vC[0]*cos2_theta[i] + (vC[1] + vC[2])*cossin_theta[i] + vC[3]*sin2_theta[i]);
        //        lambda_ens[i] = vC[0]*cos2_theta[i] + (vC[1] + vC[2])*cossin_theta[i] + vC[3]*sin2_theta[i];
        // find the weighted ensemble stress at each theta
        Sf[i] = odfweight[i] * ensemblemodel.stress(lambda_ens[i], lambda_recruit_scaling[i]);
    }
    
    {
        double sum11 = 0.0;
        double sum12 = 0.0;
        double sum22 = 0.0;
        
        for (int i = 0; i < 31; i++) {
            sum11 += Sf[i] * cos2_theta[i];
            sum12 += Sf[i] * cossin_theta[i];
            sum22 += Sf[i] * sin2_theta[i];
        }
        
        res[0] = this->kC * sum11;
        res[1] = this->kC * sum12;
        res[2] = this->kC * sum12;
        res[3] = this->kC * sum22;
    }
    
    return 0;
}







