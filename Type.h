//
//  Type.h
//  FinanceTools
//
//  Created by Alexandre HUMEAU on 28/01/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef FinanceTools_Type_h
#define FinanceTools_Type_h

#include <iostream>
#include <vector>

typedef double Time;
typedef double Rate;

typedef char* Err;

typedef std::size_t Integer;
typedef double Decimal;
typedef double Real;

typedef std::vector<double> DVector;
typedef std::vector<std::vector<double> > DMatrix;

#ifndef EPSILON
#define EPSILON 0.0000000001 // 10^(-10)
#endif

#ifndef MAX_DOUBLE
#define MAX_DOUBLE 1.0e308 //std::numeric_limits<double>::infinity
#endif

#endif
