//
//  MathFunctions.h
//  Seminaire
//
//  Created by Alexandre HUMEAU on 04/03/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef Seminaire_MathFunctions_h
#define Seminaire_MathFunctions_h

#include <cmath>
#include "TermStructure.h"
#include "Option.h"

namespace MathFunctions {
    double NormalCumulativeDistribution(const double x);

    // function to return sqrt(PI * x) --> very useful
    double sqrtpi(double x);

    // More accurate cumulative normal functions
    double AccCumNorm(double x);
    double AccBivarCumNorm(double x, double y, double rho);
    double InvCumNorm(double p);
    
    //  Beta function Ornstein Ulhenbeck process
    
#ifndef BETAOUTHRESHOLD
#define BETAOUTHRESHOLD 1.e-07
#endif
    
    double Beta_OU(double dLambda, double dt);
    
    //  Function to compute sum(exp(dLambda * u)du, u=dt1..dt2)
    double SumExp(double dLambda, double dt1, double dt2) ;
    
	    // Black-Scholes Function
    double BlackScholes(double dForward, double dStrike, double dStdDev, Finance::OptionType eOptionType);


}
#endif
