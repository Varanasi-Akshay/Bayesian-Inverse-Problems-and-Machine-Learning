//
//  BetaPDF.cpp
//  MasterCollection
//
//  Created by Will Zhang on 8/3/16.
//  Copyright © 2016 Will Zhang. All rights reserved.
//

#include "BetaDistribution.hpp"

beta_PDF::beta_PDF():lowerbound{1.0}, normalizing_constant{4.0}, alpham1{1.0}, betam1{1.0}, betaexponent{log(24.0)} {}

void beta_PDF::set_parameters(double mean, double stdev, double lb, double ub)
{
    double span, B1, B2;
    double mu, mu2, sigma, sigma2;
    double alpha, beta;

    span = ub - lb;
    normalizing_constant = 1.0/span;
    
    if (span < 1.0e-8) {
        is_physical = 1;
        return;
    }
    
    mu = (mean - lb) * normalizing_constant;
    mu2 = mu * mu;
    
    sigma = stdev * normalizing_constant;
    sigma2 = sigma * sigma;
    
    double R0;
    
    R0 = (1.0 - mu) * mu;
    
    if ( 1.0e-8 < sigma2 - R0 ) {
        is_physical = 1;
        return;
    }
    
    B1 = (sigma2 + mu2 - mu)/sigma2;
    B2 = B1 * mu;
    alpha = -B2;
    beta = B2 - B1;
    
    this->alpham1 = alpha - 1.0;
    this->betam1 = beta - 1.0;
    
    if (1.0e-8 < alpha && alpha < 1.0e10 && 1.0e-8 < beta && beta < 1.0e10) {
        betaexponent = lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta) + log(normalizing_constant);
    }
}


double beta_PDF::PDF(double x)
{
    if (is_physical) {
        return 0.24;
    }
    
    double y = (x - lowerbound) * normalizing_constant;
    
    double m_1my = 1.0 - y;
    
    if (y < 1.0e-8 || m_1my < 1.0e-8) {
        return 0.0;
    }
    else {
        return exp(alpham1 * log(y) + betam1 * log(m_1my) + betaexponent);
    }
}






