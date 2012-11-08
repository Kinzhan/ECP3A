//
//  MathFunctions.h
//  FinanceTools
//
//  Created by Alexandre HUMEAU on 04/03/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef FinanceTools_MathFunctions_h
#define FinanceTools_MathFunctions_h

#include <cmath>
#include "TermStructure.h"

namespace MathFunctions {
    double NormalCumulativeDistribution(const double x);

    // function to return sqrt(PI * x) --> very useful
    double sqrtpi(const double x);

    // More accurate cumulative normal functions
    double AccCumNorm(const double x);
    double AccBivarCumNorm(const double x, const double y, const double rho);
    double InvCumNorm(const double p);
    
    //  Beta function Ornstein Ulhenbeck process
    
#ifndef BETAOUTHRESHOLD
#define BETAOUTHRESHOLD 1.e-07
#endif
    
    double Beta_OU(const double dLambda, const double dt);
    
    //  Function to compute (Beta_OU(dLambda, dt1) - Beta_OU(dLambda, dt2)) * dLambda = \int_{t_1}^{t_2} 
    double SumExp(const double dLambda, const double dt1, const double dt2) ;
    
	// Determinist part for the HW1F model
	double DeterministPartZCBHW1F(const double dLambda, const double dt, const double dT, const Finance::TermStructure <double, double> & dTermStructure);
    
    // Black-Scholes Function
    double BlackScholes(const double dForward, const double dStrike, const double dStdDev, const std::size_t iPhi = 1);

}
#endif
