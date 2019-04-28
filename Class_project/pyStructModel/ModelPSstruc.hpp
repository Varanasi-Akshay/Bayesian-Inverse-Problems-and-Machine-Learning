//
//  ModelPSstruc.hpp
//  fullstructuralmodel
//
//  Created by Will Zhang on 6/12/18.
//  Copyright © 2018 Will Zhang. All rights reserved.
//

#ifndef ModelPSstruc_hpp
#define ModelPSstruc_hpp

#include <stdio.h>
#include "ModelMatrixWZ.hpp"
#include "ModelCollagenEXLstruc.hpp"


class model_ps_struc
{
protected:
    /*
      in this code it is assumed that we have two stage update scheme
                F_PS^{t-dt}       F_PS^{t-dt}
      B_0^{t-dt}    ==>    B_0^{t}    ==>   B_0^{t+dt}
         |                   |                   | 
         | F^{d-dt}          | F^{d}             | F^{d+dt}
         V                   V                   V
      B_1^{t-dt}           B_1^{t}          B_1^{t+dt}
    */
    double rate_PS;          // k in eq(30)
    
    double time_elapsed;     // s in eq(30)
    
    double kM;               // \eta_m in eq(16)
    double ratio;            // r in eq(16)
    double alpham1;          // \alpha in eq(16)
    double betam1;           // \beta in eq(16)
    
    double strainhistory[4]; // A(s^tilda) in eq(30)
    
    double kC;               // \eta_c
    double ODFmean;          // fibre orientation func
    double ODFstdev;
    double Dxmean;           // recuitment func
    double Dxstdev;
    double Dxlb;
    double Dxub;
    double Dx_aniso;         // anisotropy of recuitment: ortho / preffered
    double kI;               // \eta_I
    
    double F11_pushforward;  // think it as F_PS^(t-dt) if vF is F_PS^(t)
    double F12_pushforward;  // and strainhistory is A(t)
    double F21_pushforward;
    double F22_pushforward;
    
public:
    
    double fraction_old;    // exp(-k*s) in eq(30)
    double fraction_new;    // 1-exp(-k*s) to approx. the integral eq(30)
    
    model_matrix_yeoh_wz matrixmodel;
    
    // This is the updated standard structural model. 
//    model_col_tissue_struc_PF collagenmodel;
    
    // This give the full exogenously crosslinked constitutive model
    model_col_int_tissue_struc_PF collagenmodel;
    
    void calc_new_old_fraction(double rate, double time);
    
    void update_reference(double arg[4]);
    
    void update_history(double arg[4]);
    
    void updata_parameters(double arg[13]);
    
    model_ps_struc():matrixmodel(), collagenmodel(),
    kM{1.0}, ratio{1.0}, alpham1{1.0}, betam1{1.0},
    kC{100000.0}, ODFmean{0.0}, ODFstdev{0.3},
    Dxmean{1.2}, Dxstdev{0.015}, Dxlb{1.0}, Dxub{1.24}, Dx_aniso{1.0},
    kI{0.0},
    F11_pushforward{1.0}, F12_pushforward{0.0},F21_pushforward{0.0},F22_pushforward{1.0},
    strainhistory{1.0, 0.0, 0.0, 1.0}
    {calc_new_old_fraction(0.01, 0.0);};
    
    model_ps_struc(double arg[]);
    
    int stress(double vF[4], double res[4], bool crosslinked, bool withoutmatrix);
    
};


#endif /* ModelPSstruc_hpp */
