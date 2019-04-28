//
//  ModelCollagenEnsemble.cpp
//  fullstructuralmodel
//
//  Created by Will Zhang on 5/11/18.
//  Copyright © 2018 Will Zhang. All rights reserved.
//



#include "ModelCollagenEnsemble.hpp"


void model_col_ens_struc_PF::set_parameters(double mu, double sigma, double lb, double ub)
{
    mean = mu;
    stdev = sigma;
    
    lower_limit = lb;
    upper_limit = ub;
    
    Dx_PDF.set_parameters(mu, sigma, lb, ub);
}



model_col_ens_struc_PF::model_col_ens_struc_PF(): mean(1.2), stdev(0.015), lower_limit(1.0), upper_limit(1.24), Dx_PDF(1.2, 0.12, 1.12, 1.24), is_deltafunction(0)
{
};


model_col_ens_struc_PF::model_col_ens_struc_PF(double mu, double sigma, double lb, double ub):
mean(mu), stdev(sigma), lower_limit(lb), upper_limit(ub), Dx_PDF(mu, sigma, lb, ub), is_deltafunction(0)
{
}





double model_col_ens_struc_PF::stress(double lambda_fiber, double lambda_pushforward)
{
    
    double R0; //Dummy variable
    
    double ll_of_integration, ul_of_integration, span;
    
    double x1_t[gq_n], x0_t[gq_n], PDF_recruitment[gq_n];
    
    double m_1_lambdapushforward;
    
    model_col_fiber_PF fibermodel(lambda_fiber);    //fiber level model response for a give stretch
    

    /* Compute the limits of integration*/
    m_1_lambdapushforward = 1.0 / lambda_pushforward;
    
    R0 = lower_limit * m_1_lambdapushforward;
    ll_of_integration = fmax(R0, 1.0); // Calculating the lower limit (ll) for integration
    
    R0 = upper_limit * m_1_lambdapushforward;
    ul_of_integration = fmin(R0, lambda_fiber); // Calculating the upper limit (ul) for integration
    
    span = ul_of_integration - ll_of_integration;
    
    
    
    /* Compute the gauss quad integration points in beta_0 and beta_1 */
    
    for (int i = 0; i < gq_n; i++) {
        x1_t[i] = span * gq_abs[i] + ll_of_integration; // Compute the guass quadrature integration points in the beta_1 state
        x0_t[i] = x1_t[i] * lambda_pushforward; // Compute the guass quandrature integration points in the original beta_0 state
    }
    
    /* Compute the recruitment PDF function D(x) in the original state beta_0 */
    
    for (int i = 0; i < gq_n; i++) {
        PDF_recruitment[i] = Dx_PDF.PDF(x0_t[i]);
    }
    
    /* Integration over the recruitment function for the ensemble stress */
    
    double sum = 0.0;
    
    for (int i = 0; i < gq_n; i++) {
        sum = sum + gq_weight[i] * PDF_recruitment[i] * fibermodel.stress(x1_t[i]);
    }
    
    return (span * lambda_pushforward) * sum; // span is used to map the weight given the integration range, where are lambda_pushforward is used to scale the push forwarded distribution
}


double model_col_ens_struc_PF::energydensity(double lambda_fiber, double lambda_pushforward)
{
    return 0.0;     //  Currently in a dummy function will be filled in at a later time.
}
