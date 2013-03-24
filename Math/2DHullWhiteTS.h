/*
 *  2DHullWhiteTS.h
 *  Seminaire
 *
 *  Created by Emile on 1/23/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef Seminaire_2DHullWhiteTS_h
#define Seminaire_2DHullWhiteTS_h

#include "Integral.h"
#include "Require.h"

namespace Maths {
    class TwoDimHullWhiteTS : public TermStructureIntegral
    {
    protected:
        //double dLambda1_;
		//double dLambda2_;
    public:
        TwoDimHullWhiteTS(const Finance::TermStructure<double,double> & sTermStructure1, const Finance::TermStructure<double,double> & sTermStructure2/*, const double dLambda1, const double dLambda2*/);
        virtual ~TwoDimHullWhiteTS();
        
        //  primitive of function f() define in the above integral
        //  this function now computes : 
        // \int_{A}^{B} (1 - e^{-\lambda_1 (S_1 - u)}) (1 - e^{-\lambda_2 (S_2 - u)}) du / (\lambda_1 * \lambda_2) 
        // (20/02/2013 A.H. bug fix)
		virtual double TwoDimSubIntegral(double dA, double dB, double dS1, double dS2, double dLambda1, double dLambda2) const ;
		//  overloading the Integral method to compute the integrals of the sum
        //  This function now computes \int_{T1}^{T2} \Gamma_1(u,S_1) \Gamma_2(u, S_2) du (20/02/2013 A.H. bug fix)
        virtual double Integral(double dT1, double dT2, double dS1, double dS2, double dLambda1, double dLambda2) const;
		virtual double SubIntegral(double dA, double dB) const ;
    };
}

#endif
