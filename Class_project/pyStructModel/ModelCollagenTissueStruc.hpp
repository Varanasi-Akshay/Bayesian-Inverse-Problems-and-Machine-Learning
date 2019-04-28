//
//  ModelCollagenTissueStruc.hpp
//  fullstructuralmodel
//
//  Created by Will Zhang on 5/11/18.
//  Copyright © 2018 Will Zhang. All rights reserved.
//

#ifndef ModelCollagenTissueStruc_hpp
#define ModelCollagenTissueStruc_hpp



#include "BetaDistribution.hpp"
#include "ModelCollagenEnsemble.hpp"


class model_col_tissue_struc_PF
{
    //  The parameters of this model are
    //      kC = arg[0]
    //      ODFmean = arg[1]
    //      ODFstdev = arg[2]
    //      Dxmean = arg[3]
    //      Dxstdev = arg[4]
    //      Dxlb = arg[5]
    //      DXub = arg[6]
    //      Dxaniso = arg[7]
    //      kI = arg[8]
    //
    //  The constants of this model is
    //      PF11 = arg[9]
    //      PF12 = arg[10]
    //      PF21 = arg[11]
    //      PF22 = arg[12]
    //
    //  The argument of this model are
    //      F11, F12, F21, F22
    
protected:
    
    
    double F11_pushforward;
    double F12_pushforward;
    double F21_pushforward;
    double F22_pushforward;
    
    const int gq_n = 31;    // The number of gauss quadrature points useed
    const double *gq_abs = C_11_A31;    // The non transformed abscissas
    const double *gq_weight = C_11_W31; // The non-transformed weights
    
    double cos2_theta[31];
    double sin2_theta[31];
    double cossin_theta[31];
    double odfweight[31];
    
    double lambda_pushforward[31];
    double C_pushforward[31];
    
    double lambda_odf_anisotropy[31];
    
    double lambda_recruit_scaling[31];
    
    double int_range_theta; //range of integration over the splay
    
    double m_1_det;
    
    void compute_range(double stdev){int_range_theta = fmin(M_PI_2, 4.0 * stdev);}
    
public:
    
    double kC;
    double ODFmean;
    double ODFstdev;
    double Dxmean;
    double Dxstdev;
    double Dxlb;
    double Dxub;
    double kI;
    double Dx_aniso;
    
    beta_ODF fiber_odf;
    model_col_ens_struc_PF ensemblemodel;
    
    void set_parameters();
    
    model_col_tissue_struc_PF(): kC{100000.0}, kI{2000.0}, fiber_odf(0.0, 0.3),
    ensemblemodel(1.2, 0.015, 1.0, 1.24), Dx_aniso{1.0}, int_range_theta(M_PI_2), m_1_det{1.0}
    {
        double theta; //    orientation in the current reference state
        double omega; //    orientation in the original reference state
        double cos_omega;
        double sin_omega;
        double dummie_cos;
        double dummie_sin;
        
        for (int i = 0; i < 31; i++) {
            theta = int_range_theta * C_11_A31[i];
            dummie_cos = cos(theta);
            dummie_sin = sin(theta);

            cos2_theta[i] = dummie_cos*dummie_cos;
            sin2_theta[i] = dummie_sin*dummie_sin;
            cossin_theta[i] = dummie_cos*dummie_sin;
            
            omega = theta;
            cos_omega = cos(omega);
            sin_omega = sin(omega);
            
            lambda_odf_anisotropy[i] = 1.0;
            
            C_pushforward[i] = 1.0;
            lambda_pushforward[i] = 1.0;
            
            lambda_recruit_scaling[i] = 1.0;
            
            odfweight[i] = int_range_theta * C_11_W31[i] * fiber_odf.ODF(omega);
        }
    };
    
    model_col_tissue_struc_PF(double arg[13]):
    kC{arg[0]},
    ODFmean{arg[1]},
    ODFstdev{arg[2]},
    fiber_odf(arg[1], arg[2]),
    ensemblemodel(arg[3], arg[4], arg[5], arg[6]),
    Dx_aniso{arg[7]},
    kI{arg[8]},
    F11_pushforward{arg[9]},
    F12_pushforward{arg[10]},
    F21_pushforward{arg[11]},
    F22_pushforward{arg[12]},
    m_1_det{1.0/(arg[9]*arg[12] - arg[10]*arg[11])}
    {
        compute_range(arg[2]);
        double theta;
        double omega; //    orientation in the original reference state
        double cos_omega;
        double sin_omega;
        double dummie_cos;
        double dummie_sin;
        
        double mean_omega;
        
        double dummie_R0;
        double dummie_R1;
        
        mean_omega =atan2(F11_pushforward*sin(ODFmean) - F21_pushforward*cos(ODFmean),
                          F22_pushforward*cos(ODFmean) - F12_pushforward*sin(ODFmean));
        
        for (int i = 0; i < 31; i++) {
            theta = int_range_theta * C_11_A31[i] + ODFmean;
            dummie_cos = cos(theta);
            dummie_sin = sin(theta);
            
            cos2_theta[i] = dummie_cos*dummie_cos;
            sin2_theta[i] = dummie_sin*dummie_sin;
            cossin_theta[i] = dummie_cos*dummie_sin;
            
            omega = atan2(F11_pushforward*dummie_sin - F21_pushforward*dummie_cos,
                          F22_pushforward*dummie_cos - F12_pushforward*dummie_sin);
            cos_omega = cos(omega);
            sin_omega = sin(omega);
            
            dummie_R0 = cos(omega - mean_omega);
            dummie_R1 = Dx_aniso * sin(omega - mean_omega);
            lambda_odf_anisotropy[i] = sqrt(dummie_R0*dummie_R0 + dummie_R1*dummie_R1);
            
            dummie_R0 = F11_pushforward*cos_omega + F12_pushforward*sin_omega;
            dummie_R1 = F21_pushforward*cos_omega + F22_pushforward*sin_omega;
            C_pushforward[i] = dummie_R0*dummie_R0 + dummie_R1*dummie_R1;
            lambda_pushforward[i] = sqrt(C_pushforward[i]);
            
            lambda_recruit_scaling[i] = lambda_pushforward[i]/lambda_odf_anisotropy[i];
            
            odfweight[i] = int_range_theta * C_11_W31[i] * fiber_odf.ODF(omega) * C_pushforward[i] * m_1_det;
        }
    };
    
    int stress(double vF[4], double res[4]);
    
    
};


#endif /* ModelCollagenTissueStruc_hpp */
