//
//  FXMultiCurve.cpp
//  Seminaire
//
//  Created by Alexandre HUMEAU on 12/12/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "FXMultiCurve.h"
#include "MathFunctions.h"
#include "ForwardRate.h"

namespace Processes {
    
    FXMultiCurve::FXMultiCurve()
    {}
    
    FXMultiCurve::FXMultiCurve(const Finance::TermStructure<double, double> & sFXVolatilityTS, const Processes::LinearGaussianMarkov & sForeignEconomy, const Processes::LinearGaussianMarkov & sDomesticEconomy, const Finance::TermStructure<double, double> & sForeignFXCorrelation)
    {
        sFXVolatilityTS_            = sFXVolatilityTS;
        sForeignEconomy_            = sForeignEconomy;
        sDomesticEconomy_           = sDomesticEconomy;
        sForeignFXCorrelation_      = sForeignFXCorrelation;
    }
    
    FXMultiCurve::~FXMultiCurve()
    {}
    
    double FXMultiCurve::QuantoAdjustmentAdditive(const double dT1, const double dT2)
    {
        Finance::ForwardRate sForeignForwardRate(sForeignEconomy_.GetYieldCurve());
        return sForeignForwardRate.FwdRate(dT1, dT2) * (QuantoAdjustmentMultiplicative(dT1, dT2) - 1);
    }
    
    double FXMultiCurve::QuantoAdjustmentElementary(const double dStart, const double dEnd) const
    {
        double  dSigmaForeignRate = sForeignEconomy_.GetSigma().Interpolate(dStart),
                dVolFX = sFXVolatilityTS_.Interpolate(dStart),
                dRhoForeignFX = sForeignFXCorrelation_.Interpolate(dStart);
        //  We assume that between dStart and dEnd the three values of term structures are constant
        return dSigmaForeignRate * dVolFX * dRhoForeignFX * MathFunctions::SumExp(sForeignEconomy_.GetLambda(), dStart, dEnd);
    }
    
    double FXMultiCurve::QuantoAdjustmentMultiplicative(const double dT1, const double dT2)
    {
        if (sFXVolatilityTS_.IsSameTermStructure(sForeignFXCorrelation_) && sFXVolatilityTS_.IsSameTermStructure(sForeignEconomy_.GetSigma()))
        {
            // Compute \int_{t}^{dT1} a_f(u) \sigma_X(u) \rho_{f,X}(u) du
            
            //  Term structure in FX volatility, Foreign Economy and Foreign Eco & FX correlation are the same
            std::size_t iNTS = sFXVolatilityTS_.GetNbVariables();
            double dStartTS = sFXVolatilityTS_.GetVariables()[0];
            
            //  First integration between 0 and dStartTS
            double dResult = QuantoAdjustmentElementary(0, dStartTS);
            
            //  Loop over the term structure
            std::size_t iTS;
            for (iTS = 0 ; iTS < iNTS - 1 ; ++iTS)
            {
                double dStartLoc = sFXVolatilityTS_.GetVariables()[iTS], dEndLoc = sFXVolatilityTS_.GetVariables()[iTS + 1];
                if (dEndLoc <= dT1)
                {
                    dResult += QuantoAdjustmentElementary(dStartLoc, dEndLoc);
                }
                else
                {
                    break;
                }
            }
            
            
            //  Last integration between dEndTS and dT1
            double dEndTS = sFXVolatilityTS_.GetVariables()[iTS];
            dResult += QuantoAdjustmentElementary(dEndTS, dT1);
            
            //  End of calculation of integral
            double dLambdaForeign = sForeignEconomy_.GetLambda();
            return exp((MathFunctions::Beta_OU(dLambdaForeign, dT2) - MathFunctions::Beta_OU(dLambdaForeign, dT1)) * dResult);
        }
        else
        {
            //  Merge term structure
            sForeignFXCorrelation_.MergeTermStructure(sFXVolatilityTS_);
            Finance::TermStructure<double, double> sForeignEconomyVolatility = sForeignEconomy_.GetSigma();
            sForeignFXCorrelation_.MergeTermStructure(sForeignEconomyVolatility);
            sFXVolatilityTS_.MergeTermStructure(sForeignEconomyVolatility);
            sForeignEconomy_.SetSigma(sForeignEconomyVolatility);
            
            //  Call function QuantoAdjustmentMultiplicative
            return QuantoAdjustmentMultiplicative(dT1, dT2);
        }
    }
}