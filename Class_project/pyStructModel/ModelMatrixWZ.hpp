//
//  ModelMatrixWZ.hpp
//  fullstructuralmodel
//
//  Created by Will Zhang on 5/14/18.
//  Copyright © 2018 Will Zhang. All rights reserved.
//

#ifndef ModelMatrixWZ_hpp
#define ModelMatrixWZ_hpp



class model_matrix_yeoh_wz
{
private:
    double kM;
    double ratio;
    double alpham1;
    double betam1;
    
    double F11_pushforward;
    double F12_pushforward;
    double F21_pushforward;
    double F22_pushforward;
    
    double det_pushforward;
    
    double vFinv_pushforward[4];
    double vBinv_pushforward[4];
    
public:
    
    model_matrix_yeoh_wz();
    
    model_matrix_yeoh_wz(double arg[4]);
    
    void calc_tensorinverse(double vT[4], double vTinv[4], double &det);
    
    double calc_invariant_1(double vReferencestate[4], double vC[4]);
    
    int stress(double strainhistoryA[4], double vF[4], double res[4]);
    
};


#endif /* ModelMatrixWZ_hpp */
