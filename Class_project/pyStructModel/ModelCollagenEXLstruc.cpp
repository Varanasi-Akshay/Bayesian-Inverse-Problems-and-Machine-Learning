//
//  ModelCollagenEXLstruc.cpp
//  fullstructuralmodel
//
//  Created by Will Zhang on 5/15/18.
//  Copyright © 2018 Will Zhang. All rights reserved.
//


#include "kinematics.hpp"
#include "ModelCollagenEXLstruc.hpp"



int model_col_int_tissue_struc_PF::int_stress(double vF[4], double res[4]){
    
    double vC[4];
    
    double lambda_ens[31];
    
    double Sf[31][31];
    
    //  Calculate the right Cauchy Green Strain tensor
    calc_C_from_F(vF, vC);
    
    //  Calculate the ensemble stretches
    for (int i = 0; i < 31; i++) {
        lambda_ens[i] = sqrt(vC[0]*cos2_theta[i] + (vC[1] + vC[2])*cossin_theta[i] + vC[3]*sin2_theta[i]);
    }
    
    //  Calculate the ensemble stresses
    for (int i = 0; i < 31; i++) {
        for (int j = 0; j < 31; j++) {
            Sf[i][j] = odfweight[i] * odfweight[j] * ensembleinteraction.stress(lambda_ens[i], lambda_recruit_scaling[i], lambda_ens[j], lambda_recruit_scaling[j]);
        }
    }
    
    {
        double sum11 = 0.0;
        double sum12 = 0.0;
        double sum22 = 0.0;
        
        double scaling;
        
        for (int i = 0; i < 31; i++) {
            for (int j = 0; j < 31; j++) {
                scaling = lambda_ens[i]/lambda_ens[j];
                sum11 += Sf[i][j] * (cos2_theta[i] / scaling + cos2_theta[i] * scaling);
                sum12 += Sf[i][j] * (cossin_theta[i] / scaling + cossin_theta[i] * scaling);
                sum22 += Sf[i][j] * (sin2_theta[i] / scaling + sin2_theta[i] * scaling);
            }
        }
        
        res[0] = this->kC * sum11;
        res[1] = this->kC * sum12;
        res[2] = this->kC * sum12;
        res[3] = this->kC * sum22;

    }
    
    
    return 0;
}





int model_col_int_tissue_struc_PF::stress_total(double vF[4], double res[4])
{
    double collagenstress[4];
    double interactionstress[4];
    
    int errcol, errint;
    
    errcol = stress(vF, collagenstress);
    errint = int_stress(vF, interactionstress);
    
    for (int i = 0; i<4; i++) {
        res[i] = collagenstress[i] + interactionstress[i];
    }
    
    return 0;
}












