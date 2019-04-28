//
//  main.cpp
//  testmain
//
//  Created by Will Zhang on 11/29/17.
//  Copyright © 2017 Will Zhang. All rights reserved.
//

// #include <iostream>
#include "stress.hpp"


void stress(double parameters[23], double strain[4], double res[4], bool crosslinked, bool withoutmatrix) {
    
    // double parameters[23] = {100, 0.0, 1.0 , 1.0, 100000, 0.0, 0.3, 1.2, 0.015, 1.0, 1.24, 1.0, 1000, 1.0, 0.0, 0.0, 1.0, 0.0001, 1000, 1.05, 0.0, 0.0, 1.05};

    model_ps_struc model(parameters);

    // double strain[4] = {1.2, 0.0, 0.0, 1.4};
    // double res[4];

//    res = model.strainenergy(strain);

  
    int err; 
    err = model.stress(strain, res, crosslinked, withoutmatrix);

    return;
}
