//
//  CollagenTissueModels.cpp
//  fullstructuralmodel
//
//  Created by Will Zhang on 11/29/17.
//  Copyright © 2017 Will Zhang. All rights reserved.
//



#include "ModelCollagenTissueExp.hpp"

int col_tissue_PF_exp::stress(double vC[4], double res[4])
{
    
    double lambda_ens[31];
    double Sf[31]; //Ensemble stress at each theta
    
    for (int i = 0; i < 31; i++) {
        lambda_ens[i] = sqrt(vC[0]*cos2_theta[i] + (vC[1] + vC[2])*cossin_theta[i] + vC[3]*sin2_theta[i]);
//        lambda_ens[i] = vC[0]*cos2_theta[i] + (vC[1] + vC[2])*cossin_theta[i] + vC[3]*sin2_theta[i];
        // find the weighted ensemble stress at each theta
        Sf[i] = odfweight[i] * ensmodel.stress(lambda_ens[i]);
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


double col_tissue_PF_exp::strainenergy(double vC[3])
{
    double lambda_ens[31];
    double Psif[31]; //Ensemble stress at each theta
    
    for (int i = 0; i < 31; i++) {
        lambda_ens[i] = sqrt(vC[0]*cos2_theta[i] + (vC[1] + vC[2])*cossin_theta[i] + vC[3]*sin2_theta[i]);
//        lambda_ens[i] = vC[0]*cos2_theta[i] + (vC[1] + vC[2])*cossin_theta[i] + vC[3]*sin2_theta[i];
        // find the weighted strain energy at each theta
        Psif[i] = odfweight[i] * ensmodel.strainenergy(lambda_ens[i]);
    }
    
    double sumer = 0.0;
    
    for (int i = 0; i < 31; i++) {
        sumer += Psif[i];
    }
    
    return this->kC * sumer;
}
