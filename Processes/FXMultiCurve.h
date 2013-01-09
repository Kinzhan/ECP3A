//
//  FXMultiCurve.h
//  Seminaire
//
//  Created by Alexandre HUMEAU on 12/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef Seminaire_FXMultiCurve_h
#define Seminaire_FXMultiCurve_h

#include "TermStructure.h"
#include "HullWhite.h"

namespace Processes {
    //  Definition of a class with an FX Analogy in order to have the multi-curve framework with the stochastic basis spread (see Bianchetti 2009 (Two Curves, one price : pricing and hedging interest rate derivatives using different yield curves for discounting and forwarding
    class FXMultiCurve
    {
    public:
        FXMultiCurve();
        FXMultiCurve(const Finance::TermStructure<double, double> & sFXVolatilityTS, const Processes::LinearGaussianMarkov & sForeignEconomy, const Processes::LinearGaussianMarkov & sDomesticEconomy, const Finance::TermStructure<double, double> & sForeignFXCorrelation);
        virtual ~FXMultiCurve();
        
        // Computation of the integral \int_{dStart}^{dEnd} a(s) \sigma_X(u) \rho_{f,X}(u) du where all termstructure are constant between dStart and dEnd
        virtual double QuantoAdjustmentElementary(const double dStart, const double dEnd) const;
        
        virtual double QuantoAdjustmentMultiplicative(const double dT1, const double dT2);
        virtual double QuantoAdjustmentAdditive(const double dStart, const double dEnd);
        
    protected:
        Finance::TermStructure<double, double> sFXVolatilityTS_;
        Processes::LinearGaussianMarkov sForeignEconomy_, //     Forwarding curve
                                        sDomesticEconomy_; //    Discount Curve
        Finance::TermStructure<double, double> sForeignFXCorrelation_;
    };
}

#endif
