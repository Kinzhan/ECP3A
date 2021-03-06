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
#include "HullWhiteTSCorrection.h"
#include "HullWhiteTS.h"
#include "Weights.h"

namespace Processes {
    StochasticBasisSpread::StochasticBasisSpread()
    {}
    
    StochasticBasisSpread::~StochasticBasisSpread()
    {}
    
    double StochasticBasisSpread::CorrelationSpreadOIS(const Finance::TermStructure<double, double> & sSigmaOISTS, 
                                                       const Finance::TermStructure<double, double> & sSigmaCollatTS, 
                                                       double dLambdaOIS, 
                                                       double dLambdaCollat, 
                                                       double dRhoCollatOIS,
                                                       double dt,
                                                       double dT) const
    {
        return (dRhoCollatOIS * sSigmaCollatTS.Interpolate(dt) * exp(-dLambdaCollat * (dT - dt)) - sSigmaOISTS.Interpolate(dt) * exp(-dLambdaOIS * (dT - dt))) / VolSpread(sSigmaOISTS, sSigmaCollatTS, dLambdaOIS, dLambdaCollat, dRhoCollatOIS, dt, dT);
    }
    
    double StochasticBasisSpread::VolSpread(const Finance::TermStructure<double, double> & sSigmaOISTS, 
                                            const Finance::TermStructure<double, double> & sSigmaCollatTS, 
                                            double dLambdaOIS, 
                                            double dLambdaCollat, 
                                            double dRhoCollatOIS,
                                            double dt,
                                            double dT) const
    {
        Utilities::require(dRhoCollatOIS >= -1.0 && dRhoCollatOIS <= 1, "Correlation between Collat and OIS is not between -1 and 1");
        double dVolCollat = sSigmaCollatTS.Interpolate(dt), dVolOIS = sSigmaOISTS.Interpolate(dt);
        double dVariance = dVolCollat * dVolCollat * exp(-2.0 * dLambdaCollat * (dT - dt)) + dVolOIS * dVolOIS * exp(-2.0 * dLambdaOIS * (dT - dt)) - 2 * dRhoCollatOIS * dVolOIS * dVolCollat * exp(-(dLambdaOIS + dLambdaCollat) * (dT - dt));
        return sqrt(dVariance);
    }
    
    double StochasticBasisSpread::LiborQuantoAdjustmentMultiplicative(Finance::TermStructure<double, double> & sSigmaOISTS, 
                                                                 Finance::TermStructure<double, double> & sSigmaCollatTS, 
                                                                 double dLambdaOIS, 
                                                                 double dLambdaCollat, 
                                                                 double dRhoCollatOIS,
                                                                 double dt,
                                                                 double dT1,
                                                                 double dT2,
                                                                 std::size_t iNIntervals) const
    {
        //  numeric integration for now (may need exact computation)
        double dResult = 0;
        
		Maths::HullWhiteTS sCollatHWTS(sSigmaCollatTS, dLambdaCollat), sOISHWTS(sSigmaOISTS, dLambdaOIS);
        for (std::size_t iInterval = 0 ; iInterval < iNIntervals ; ++iInterval)
        {
            //  Middle of each small interval
            double dMiddle = dt + (iInterval + 0.5) * (dT1 - dt) / iNIntervals;
            double df = (sCollatHWTS.Integral(dMiddle, dT2) - sCollatHWTS.Integral(dMiddle, dT1)) * (dRhoCollatOIS * sOISHWTS.Integral(dMiddle, dT2) - sCollatHWTS.Integral(dMiddle, dT2));
            dResult += df * (dT1 - dt);
        }
        dResult /= iNIntervals;
        return exp(-dResult);
    }
    
    double StochasticBasisSpread::SwapQuantoAdjustmentMultiplicative( const Finance::TermStructure<double, double> & sSigmaOISTS, 
																	  const Finance::TermStructure<double, double> & sSigmaCollatTS,
                                                                      double dLambdaOIS,
                                                                      double dLambdaCollat, 
                                                                      double dRhoCollatOIS,
																	  const Finance::YieldCurve & sYieldCurveOIS,
																	  const Finance::YieldCurve & sYieldCurveCollat,
                                                                      double dt,
																	  const std::vector <double> dS,
																	  const std::vector <double> dT) const
    {
		Finance::Weights sWeightsOIS(sYieldCurveOIS, dS) ;
		Finance::Weights sWeightsCollat(sYieldCurveCollat, dS) ;
		std::vector <double> dWeightsOIS = sWeightsOIS.GetWeights() ;
		std::vector <double> dWeightsCollat = sWeightsCollat.GetWeights() ;
		
		std::size_t iSizeT = dT.size() ;
		std::size_t iSizeS = dS.size() ;
		
		double dT_0 = dT[0], dT_n = dT[iSizeT-1], dResult = 0.0 ;
		Finance::DF sDFCollat(sYieldCurveCollat) ;
		double dDF_0 = sDFCollat.DiscountFactor(dT_0) ;
		double dDF_n = sDFCollat.DiscountFactor(dT_n) ;
		double dDFratio_1 = dDF_0 / (dDF_0 - dDF_n) ;
		double dDFratio_2 = dDF_n / (dDF_0 - dDF_n) ;
		
		Maths::TwoDimHullWhiteTS sTwoDimHullWhiteTSdf(sSigmaOISTS, sSigmaCollatTS), sTwoDimHullWhiteTSff(sSigmaCollatTS, sSigmaCollatTS) ;
		
		// first term, cf. report
		for (std::size_t iOIS = 0; iOIS < iSizeS - 1; ++iOIS) {
			for (std::size_t iCollat = 0; iCollat < iSizeS - 1; ++iCollat) {
				dResult -=  dRhoCollatOIS * dWeightsOIS[iOIS] * dWeightsCollat[iCollat] * sTwoDimHullWhiteTSdf.Integral(dt, dT_0, dS[iOIS + 1], dS[iCollat + 1], dLambdaOIS, dLambdaCollat) ;
			}
		}
		
		// second term
		for (std::size_t iCollat_1 = 0; iCollat_1 < iSizeS - 1; ++iCollat_1) {
			for (std::size_t iCollat_2 = 0; iCollat_2 < iSizeS - 1; ++iCollat_2) {
				dResult +=  dWeightsCollat[iCollat_1] * dWeightsCollat[iCollat_2] * sTwoDimHullWhiteTSff.Integral(dt, dT_0, dS[iCollat_1 + 1], dS[iCollat_2 + 1], dLambdaCollat, dLambdaCollat) ;
			}
		}
        
        // third term
		for (std::size_t iOIS = 0; iOIS < iSizeS - 1; ++iOIS) {
			dResult +=  dRhoCollatOIS * dWeightsOIS[iOIS] * dDFratio_1 * sTwoDimHullWhiteTSdf.Integral(dt, dT_0, dT_0, dS[iOIS + 1], dLambdaCollat, dLambdaOIS) ; // to check
		}
		
		// fourth term
		for (std::size_t iCollat = 0; iCollat < iSizeS - 1; ++iCollat) {
			dResult -=  dWeightsCollat[iCollat] * dDFratio_1 * sTwoDimHullWhiteTSff.Integral(dt, dT_0, dT_0, dS[iCollat + 1], dLambdaCollat, dLambdaCollat) ;
		}
        
        // fifth term
		for (std::size_t iOIS = 0; iOIS < iSizeS - 1; ++iOIS) {
			dResult -=  dRhoCollatOIS * dWeightsOIS[iOIS] * dDFratio_2 * sTwoDimHullWhiteTSdf.Integral(dt, dT_0, dT_n, dS[iOIS + 1], dLambdaCollat, dLambdaOIS) ;
		}
		
		// sixth term
		for (std::size_t iCollat = 0; iCollat < iSizeS - 1; ++iCollat) {
			dResult +=  dWeightsCollat[iCollat] * dDFratio_2 * sTwoDimHullWhiteTSff.Integral(dt, dT_0, dT_n, dS[iCollat + 1], dLambdaCollat, dLambdaCollat) ;
		}
													
        return exp(dResult);
    }
}