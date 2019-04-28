//
//  main.cpp
//  testmain
//
//  Created by Will Zhang on 11/29/17.
//  Copyright © 2017 Will Zhang. All rights reserved.
//

#ifndef STRESSLIB_H
#define STRESSLIB_H

// #include <iostream>
#include "BetaDistribution.hpp"
#include "ModelCollagenTissueExp.hpp"
#include "ModelPSstruc.hpp"

void stress(double parameters[23], double strain[4], double res[4], bool crosslinked, bool withoutmatrix); 

#endif
