//
//  ModelCollagenEnsembleEXLstruc.hpp
//  fullstructuralmodel
//
//  Created by Will Zhang on 5/15/18.
//  Copyright © 2018 Will Zhang. All rights reserved.
//

#ifndef ModelCollagenEnsembleEXLstruc_hpp
#define ModelCollagenEnsembleEXLstruc_hpp

#include "Gauss_Quad_Points.hpp"
#include "BetaDistribution.hpp"

class model_col_exl_ens_struc_integrand
{
private:
    
//    double lambda_alpha;
//    double lambda_beta;
    double lambda_alphaxbeta;

public:
    
    model_col_exl_ens_struc_integrand(): lambda_alphaxbeta{1.0} {};
    
    model_col_exl_ens_struc_integrand(double l_a, double l_b): lambda_alphaxbeta{l_a*l_b} {};
    
    double integrand(double x_alpha, double x_beta)
    {
        double m_1_x_alphaxbeta;
        m_1_x_alphaxbeta = 1.0/(x_alpha * x_beta);

        return m_1_x_alphaxbeta * (lambda_alphaxbeta*m_1_x_alphaxbeta - 1.0);
    };
};







class model_col_exl_ens_struc
{
    
private:
    double c0; // This is the modulus of the interaction term of the model
    
    double Dx_mean;  //  Mean of the standard deviation of the recruitment distribution
    double Dx_stdev;
    double Dx_lowerbound;  //  lower limit of integration
    double Dx_upperbound;  //  upper limit of integration
    
    beta_PDF Dx_PDF;
    
    const int gaussn = 31;
    //set quadrature abscissas
    const double *a = C_01_A31;
    const double *w = C_01_W31;
    
public:
    

    model_col_exl_ens_struc():Dx_mean{1.2}, Dx_stdev{0.015}, Dx_lowerbound{1.0}, Dx_upperbound{1.24}, Dx_PDF(1.2, 0.015, 1.0, 1.24) {};
    
    model_col_exl_ens_struc(double mu, double std, double lb, double ub):
    Dx_mean{mu}, Dx_stdev{std}, Dx_lowerbound{lb}, Dx_upperbound{ub}, Dx_PDF(mu, std, lb, ub)
    {};
    
    
    double stress(double lambda_alpha, double lambda_pf_alpha, double lambda_beta, double lambda_pf_beta);
};




#endif /* ModelCollagenEnsembleEXLstruc_hpp */





