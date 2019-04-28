//
//  ModelCollagenEXLstruc.hpp
//  fullstructuralmodel
//
//  Created by Will Zhang on 5/15/18.
//  Copyright © 2018 Will Zhang. All rights reserved.
//

#ifndef ModelCollagenEXLstruc_hpp
#define ModelCollagenEXLstruc_hpp

#include <math.h>
#include "Gauss_Quad_Points.hpp"
#include "BetaDistribution.hpp"
#include "ModelCollagenTissueStruc.hpp"
#include "ModelCollagenEnsembleEXLstruc.hpp"


class model_col_int_tissue_struc_PF: public model_col_tissue_struc_PF {

private:

    
public:
    model_col_exl_ens_struc ensembleinteraction;
    
    model_col_int_tissue_struc_PF(): model_col_tissue_struc_PF(), ensembleinteraction()  {}
    
    model_col_int_tissue_struc_PF(double arg[13]): model_col_tissue_struc_PF(arg), ensembleinteraction(arg[3], arg[4], arg[5], arg[6]) {}
    
    int int_stress(double vF[4], double res[4]);
    
    int stress_total(double vF[4], double res[4]);
    
    
};



#endif /* ModelCollagenEXLstruc_hpp */
