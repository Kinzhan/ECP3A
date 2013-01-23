//
//  StochasticBasisSpread.cpp
//  Seminaire
//
//  Created by Alexandre HUMEAU on 09/01/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "StochasticBasisSpread.h"
#include <cmath>
#include "Require.h"
#include "HullWhiteTS.h"
#include "Weights.h"

namespace Processes {
    StochasticBasisSpread::StochasticBasisSpread()
    {}
    
    StochasticBasisSpread::~StochasticBasisSpread()
    {}
    
    /*double StochasticBasisSpread::Sigma_B(double t, double T) const
    {
        //return dSigma_ * exp(-dLambda_ * (T - t));
        return 0;
    }*/
    
    double StochasticBasisSpread::CorrelationSpreadOIS(const Finance::TermStructure<double, double> & sSigmaOISTS, 
                                                       const Finance::TermStructure<double, double> & sSigmaCollatTS, 
                                                       const double dLambdaOIS, 
                                                       const double dLambdaCollat, 
                                                       const double dRhoCollatOIS,
                                                       const double dt,
                                                       const double dT) const
    {
        return (dRhoCollatOIS * sSigmaCollatTS.Interpolate(dt) * exp(-dLambdaCollat * (dT - dt)) - sSigmaOISTS.Interpolate(dt) * exp(-dLambdaOIS * (dT - dt))) / VolSpread(sSigmaOISTS, sSigmaCollatTS, dLambdaOIS, dLambdaCollat, dRhoCollatOIS, dt, dT);
    }
    
    double StochasticBasisSpread::VolSpread(const Finance::TermStructure<double, double> & sSigmaOISTS, 
                                            const Finance::TermStructure<double, double> & sSigmaCollatTS, 
                                            const double dLambdaOIS, 
                                            const double dLambdaCollat, 
                                            const double dRhoCollatOIS,
                                            const double dt,
                                            const double dT) const
    {
        Utilities::require(dRhoCollatOIS >= -1.0 && dRhoCollatOIS <= 1, "Correlation between Collat and OIS is not between -1 and 1");
        double dVolCollat = sSigmaCollatTS.Interpolate(dt), dVolOIS = sSigmaOISTS.Interpolate(dt);
        double dVariance = dVolCollat * dVolCollat * exp(-2.0 * dLambdaOIS * (dT - dt)) + dVolOIS * dVolOIS * exp(-2.0 * dLambdaCollat * (dT - dt)) - 2 * dRhoCollatOIS * dVolOIS * dVolCollat * exp(-(dLambdaOIS + dLambdaCollat) * (dT - dt));
        return sqrt(dVariance);
    }
    
    double StochasticBasisSpread::LiborQuantoAdjustmentMultiplicative(const Finance::TermStructure<double, double> & sSigmaOISTS, 
                                                                 const Finance::TermStructure<double, double> & sSigmaCollatTS, 
                                                                 const double dLambdaOIS, 
                                                                 const double dLambdaCollat, 
                                                                 const double dRhoCollatOIS,
                                                                 const double dt,
                                                                 const double dT1,
                                                                 const double dT2,
                                                                 const std::size_t iNIntervals) const
    {
        //  numeric integration for now (may need exact computation)
        double dResult = 0;
        
        Maths::HullWhiteTS sCollatHWTS(sSigmaCollatTS, dLambdaCollat), sOISHWTS(sSigmaOISTS, dLambdaOIS);
        for (std::size_t iInterval = 0 ; iInterval < iNIntervals ; ++iInterval)
        {
            //  Middle of each small interval
            double dMiddle = dt + (iInterval + 0.5) * (dT1 - dt) / iNIntervals;
            double df = exp(dLambdaCollat * dMiddle) * (exp(dLambdaOIS * dMiddle) * dRhoCollatOIS * sCollatHWTS.Integral(dMiddle, dT2) - sOISHWTS.Integral(dMiddle, dT2) * exp(dLambdaCollat * dMiddle));
            dResult += df * (dT1 - dt) / iNIntervals;
        }
        return exp(-dResult * sCollatHWTS.SubIntegral(dT1, dT2));
    }
    
    double StochasticBasisSpread::SwapQuantoAdjustmentMultiplicative(const Finance::TermStructure<double, double> & sSigmaOISTS, 
																	  const Finance::TermStructure<double, double> & sSigmaCollatTS, 
																	  const double dLambdaOIS, 
																	  const double dLambdaCollat, 
																	  const double dRhoCollatOIS,
																	  const double dt,
																	  const std::vector <double> dS,
																	  const std::vector <double> dT,
																	  const std::size_t iNIntervals) const
    {
		/*Finance::Weights sWeightsOIS(sSigmaOISTS, dS) ;
		Finance::Weights sWeightsCollat(sSigmaCollatTS, dS) ;
		
		std::size_t iSizeT = dT.size() ;
		std::size_t iSizeS = dS.size() ;
		
		double dT_0 = dT[0], dT_n = dT[iSizeT-1] ;
		
		Maths::TwoDimHullWhiteTS sTwoDimHullWhiteTS(*/
        return 0.0;
    }
}