//
//  ModelCollagenTissueExp.hpp
//  fullstructuralmodel
//
//  Created by Will Zhang on 5/11/18.
//  Copyright © 2018 Will Zhang. All rights reserved.
//

#ifndef ModelCollagenTissueExp_hpp
#define ModelCollagenTissueExp_hpp

#include "math.h"
#include "Gauss_Quad_Points.hpp"
#include "BetaDistribution.hpp"

/*****          Collagen fiber ensemble models           *****/

class col_ens_PF_exp
{
private:
    double parB;
    double lmax;
    
    double scaling;
    
public:
    
    col_ens_PF_exp(): parB(1.0), lmax(1.0), scaling{exp(- parB * (lmax - 1.0))} {}
    
    col_ens_PF_exp(double b, double max): parB{b}, lmax{max}, scaling{exp(- parB * (lmax - 1.0))} {}
    
    void set_AB(double b, double max){ this->parB = b; this->lmax = max; this->scaling = exp(- parB * (lmax - 1.0)); }
    
    double stress(double lf);
    
    double strainenergy(double lf);
};






/*****          Collagen tissue models           *****/

class col_tissue_PF_exp
{
    // The parameters are:
    //   kC     The modulus
    //   mu     The mean of the ODF
    //   sigma  The standard deviation of the ODF
    //   b      The exponent of the ensemble model
    //   lfmax  The maximum ensemble stretch applyed, this is for scaling to improve the covariance between parameters This value is set to 1.0 normally, aka no scaling.
private:
    
    double cos2_theta[31];
    double sin2_theta[31];
    double cossin_theta[31];
    double odfweight[31];
    
    double int_range_theta; //range of integration over the splay
    
    void compute_range(double stdev){int_range_theta = fmin(M_PI_2, 4.0 * stdev);}
    
    beta_ODF fiber_odf;
    
    col_ens_PF_exp ensmodel;
    
    double kC;
    
public:
    // Default constructor
    col_tissue_PF_exp(): int_range_theta(M_PI_2), fiber_odf(0, M_PI/4.0), kC(1.0), ensmodel(1.0,1.0){
        {
            double theta;
            double hereonly_cos;
            double hereonly_sin;
            for (int i = 0; i < 31; i++) {
                theta = int_range_theta * C_11_A31[i];
                hereonly_cos = cos(theta);
                hereonly_sin = sin(theta);
                cos2_theta[i] = hereonly_cos*hereonly_cos;
                sin2_theta[i] = hereonly_sin*hereonly_sin;
                cossin_theta[i] = hereonly_cos*hereonly_sin;
                odfweight[i] = int_range_theta * C_11_W31[i] * fiber_odf.ODF(theta);
            }
        }
    }
    
    // Constructor with model parameters
    col_tissue_PF_exp(double k, double mu, double sigma, double b, double lfmax): fiber_odf(mu, sigma), kC(k),
    ensmodel(b,lfmax)
    {
        compute_range(sigma);
        
        {
            double theta;
            double hereonly_cos;
            double hereonly_sin;
            for (int i = 0; i < 31; i++) {
                theta = int_range_theta * C_11_A31[i] + mu;
                hereonly_cos = cos(theta);
                hereonly_sin = sin(theta);
                cos2_theta[i] = hereonly_cos*hereonly_cos;
                sin2_theta[i] = hereonly_sin*hereonly_sin;
                cossin_theta[i] = hereonly_cos*hereonly_sin;
                odfweight[i] = int_range_theta * C_11_W31[i] * fiber_odf.ODF(theta);
            }
        }
    }
    
    int stress(double vC[4], double res[4]);
    
    double strainenergy(double vC[4]);
};


#endif /* ModelCollagenTissueExp_hpp */
