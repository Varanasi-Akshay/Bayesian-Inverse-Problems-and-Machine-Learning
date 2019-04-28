//
//  Gauss_Quad_Points.h
//  EXL_Model
//
//  Created by Will Zhang on 8/5/15.
//  Copyright (c) 2015 Will Zhang. All rights reserved.
//
//
//#ifndef __EXL_Model__Gauss_Quad_Points__
//#define __EXL_Model__Gauss_Quad_Points__
//
//extern const double C_01_A21[];
//extern const double C_01_W21[];
//
//extern const double C_11_A31[];
//extern const double C_11_W31[];
//
//
//#endif /* defined(__EXL_Model__Gauss_Quad_Points__) */

#ifndef __EXL_Model__Gauss_Quad_Points__

#define __EXL_Model__Gauss_Quad_Points__

extern "C" {
    #ifdef __EXL_Model__Gauss_Quad_Points__
    extern const double C_01_A21[];
    extern const double C_01_W21[];
    
    extern const double C_01_A31[];
    extern const double C_01_W31[];
    
    extern const double C_11_A31[];
    extern const double C_11_W31[];

    #endif
}


#endif
