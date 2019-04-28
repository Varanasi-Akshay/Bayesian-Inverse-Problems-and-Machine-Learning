//
//  BetaODF.cpp
//  PS_Model
//
//  Created by Will Zhang on 7/29/16.
//  Copyright © 2016 Will Zhang. All rights reserved.
//

#include "BetaDistribution.hpp"



void beta_ODF::set_baseline(double stdev)
{
    double R0; //dummie variable
    
    const double wrappednormalscaling = 3.891848834996066;
    
    if (stdev > M_PI_SQRT12)
    {
        is_isotropic = 1;
    }
    else if (stdev < 1.0e-4) {
        is_isotropic = 1;
        const_value = 0.0;
    }
    else if (stdev < 0.4) {
        d = 0.0;
        const_value = 0.0;
        sigma = stdev;
    }
    else
    {
        has_baseline = 1;
        R0 = (stdev - 0.4);
        d = wrappednormalscaling * R0 * R0;
        const_value = d*M_1_PI;
        sigma = find_sigma(stdev, d);
    }
}


double beta_ODF::find_sigma(double stdev, double d)
{
    return (stdev - d * M_PI_SQRT12) / (1.0 - d);
}





double beta_ODF::find_alpha(double std)
{
    double sigma = std * M_1_PI; /* map to (0,1) */
    double sigma2 = sigma * sigma; /*saves 1 muliplication*/
    return 0.5 * (0.25 / sigma2 - 1.0); /*the distribution is always symmetric*/
}





void beta_ODF::set_betaexponent(double alpha)
{
    betaexponent = lgamma(2.0 * alpha) - 2.0 * lgamma(alpha) - M_LOG_PI;
}





double beta_ODF::ODF(double theta)
{
    double y, R0;
    
    if (is_isotropic) {
        return const_value;
    }
    
    y = (theta - mu) / M_PI + 0.5;
    R0 = y - floor(y); //dummie variable
    y = R0*(1.0 - R0);
    
    if (y < 1.0e-12) {
        return 0.0;
    }
    
    R0 = exp((alpha - 1.0) * log(y) + betaexponent);
    
    if (has_baseline) {
        return (1.0 - d) * R0 + const_value;
    }
    
    return R0;

}









