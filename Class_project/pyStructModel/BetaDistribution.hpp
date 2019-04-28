//
//  BetaPDF.hpp
//  MasterCollection
//
//  Created by Will Zhang on 8/3/16.
//  Copyright © 2016 Will Zhang. All rights reserved.
//

#ifndef BetaPDF_hpp
#define BetaPDF_hpp

#include <math.h>

class beta_PDF
{
    
private:
    
    int is_physical = 0; //bool if there is baseline, aka d > 0
    

    
public:
    
    double lowerbound;
    
    double normalizing_constant;
    
    double alpham1;
    double betam1;
    
    double betaexponent;
    
    void set_parameters(double mean, double stdev, double lb, double ub);
    
    beta_PDF();
    
    beta_PDF(double mean, double stdev, double lb, double ub):lowerbound{lb}
    { set_parameters(mean, stdev, lb, ub);}
    
    double PDF(double x);

};






const double M_PI_SQRT12 = M_PI_2/sqrt(3.0); /* Pi/sqrt(12)      */
const double M_BETAMax = 3.20484408857591466767602498042342210;   /* Beta(23/8)      */
const double M_LOG_PI = log(M_PI);   /* Beta(23/8)      */


class beta_ODF
{
    
public:
    beta_ODF(){}
    
    beta_ODF(double mean, double stdev)
    {
        mu = mean;
        set_baseline(stdev);
        alpha = find_alpha(sigma);
        alpha_m_1 = alpha - 1.0;
        set_betaexponent(alpha);
    }
    
    void set_parameters(double mean, double stdev)
    {
        mu = mean;
        set_baseline(stdev);
    }
    
    void set_baseline(double stdev);
    
    double find_sigma(double stdev, double d);
    
    double find_alpha(double std);
    
    void set_betaexponent(double alpha);
    
// Output function
    
    double ODF(double theta);
    
//private:
    
    int is_isotropic = 0; //bool if the standard deviation exceed the maximum possible PI/sqrt(12.0), or is lower than the minimum 1.0e-4 (i.e. is approximately a delta function)
    int has_baseline = 0; //bool if there is baseline, aka d > 0
    
    double const_value = M_1_PI; //returns this values if is_isotropic is true
    
    double mu = 0.0;
    double sigma = 0.3;
    double d = 0.0;
    double alpha = 13.207783890401885;
    double alpha_m_1 = 12.207783890401885;
    double betaexponent = 17.180575633615234;
    
    
    
};



#endif /* BetaPDF_hpp */
