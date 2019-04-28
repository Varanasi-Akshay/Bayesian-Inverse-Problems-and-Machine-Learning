//
//  ModelCollagenEnsembleEXLstruc.cpp
//  fullstructuralmodel
//
//  Created by Will Zhang on 5/15/18.
//  Copyright © 2018 Will Zhang. All rights reserved.
//

#include "ModelCollagenEnsembleEXLstruc.hpp"




double model_col_exl_ens_struc::stress(double lambda_alpha, double lambda_pf_alpha, double lambda_beta, double lambda_pf_beta)
{
    
    double R0; //Dummy variable
    
    double ll_alpha, ul_alpha, span_alpha; // range of integration for the alpha family
    
    double ll_beta, ul_beta, span_beta; // range of integration for the beta family
    
    double x1_alpha[gaussn], x0_alpha[gaussn], PDF_recruitment_alpha[gaussn];
    
    double x1_beta[gaussn], x0_beta[gaussn], PDF_recruitment_beta[gaussn];
    
    double m_1_lambda_pf_alpha, m_1_lambda_pf_beta;
    
    model_col_exl_ens_struc_integrand integrand(lambda_alpha, lambda_beta);    //fiber level model response for a give stretch
    
    
    /* Compute the limits of integration*/
    m_1_lambda_pf_alpha = 1.0 / lambda_pf_alpha;
    m_1_lambda_pf_beta = 1.0 / lambda_pf_beta;
    
    
    /* Compute the range of integration for the alpha fiber family*/
    R0 = Dx_lowerbound * m_1_lambda_pf_alpha;
    ll_alpha = fmax(R0, 1.0); // Calculating the lower limit (ll) for integration
    R0 = Dx_upperbound * m_1_lambda_pf_alpha;
    ul_alpha = fmin(R0, lambda_alpha); // Calculating the upper limit (ul) for integration
    span_alpha = ul_alpha - ll_alpha;
    
    
    /* Compute the range of integration for the beta fiber family*/
    R0 = Dx_lowerbound * m_1_lambda_pf_beta;
    ll_beta = fmax(R0, 1.0); // Calculating the lower limit (ll) for integration
    R0 = Dx_upperbound * m_1_lambda_pf_beta;
    ul_beta = fmin(R0, lambda_beta); // Calculating the upper limit (ul) for integration
    span_beta = ul_beta - ll_beta;
    
    
    
    /* Compute the gauss quad integration points in omega_0 and omega_1 */
    
    for (int i = 0; i < gaussn; i++) {
        x1_alpha[i] = span_alpha * a[i] + ll_alpha; // Compute the guass quadrature integration points in the omega_1 state
        x0_alpha[i] = x1_alpha[i] * lambda_pf_alpha; // Compute the guass quandrature integration points in the original omega_0 state
    }
    
    for (int i = 0; i < gaussn; i++) {
        x1_beta[i] = span_beta * a[i] + ll_beta; // Compute the guass quadrature integration points in the omega_1 state
        x0_beta[i] = x1_beta[i] * lambda_pf_beta; // Compute the guass quandrature integration points in the original omega_0 state
    }
    
    /* Compute the recruitment PDF function D(x) in the original state beta_0 */
    
    for (int i = 0; i < gaussn; i++) {
        PDF_recruitment_alpha[i] = Dx_PDF.PDF(x0_alpha[i]);
        PDF_recruitment_beta[i] = Dx_PDF.PDF(x0_beta[i]);
    }
    
    /* Integration over the recruitment function for the ensemble stress */
    
    double sum = 0.0;
    
    for (int i = 0; i < gaussn; i++) {
        for (int j = 0; j < gaussn; j++) {
            sum = sum + w[i] * PDF_recruitment_alpha[i] * w[j] * PDF_recruitment_beta[j]* integrand.integrand(x1_alpha[i], x1_beta[j]);
        }
    }
    
    return (span_alpha * lambda_pf_alpha * span_beta * lambda_pf_beta) * sum; // span is used to map the weight given the integration range, where are lambda_pushforward is used to scale the push forwarded distribution
}
