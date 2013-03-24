//
//  StochasticBasisSpread.h
//  Seminaire
//
//  Created by Alexandre HUMEAU on 09/01/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#ifndef Seminaire_StochasticBasisSpread_h
#define Seminaire_StochasticBasisSpread_h

#include "TermStructure.h"
#include "2DHullWhiteTS.h"
#include "YieldCurve.h"

namespace Processes {
    //  This class model a stochastic basis spread
    class StochasticBasisSpread
    {
    protected:
    public:
        StochasticBasisSpread();
        virtual ~StochasticBasisSpread();
        
        //virtual double Sigma_B(double t, double T) const;
        virtual double CorrelationSpreadOIS(const Finance::TermStructure<double, double> & sSigmaOISTS, 
                                                           const Finance::TermStructure<double, double> & sSigmaCollatTS, 
                                                           const double dLambdaOIS, 
                                                           const double dLambdaCollat, 
                                                           const double dRhoCollatOIS,
                                                           const double dt,
                                                           const double dT) const;
        virtual double VolSpread(const Finance::TermStructure<double, double> & sSigmaOISTS, 
                                 const Finance::TermStructure<double, double> & sSigmaCollatTS, 
                                 const double dLambdaOIS, 
                                 const double dLambdaCollat, 
                                 const
								 double dRhoCollatOIS,
                                 const double dt,
                                 const double dT) const;
        
        virtual double LiborQuantoAdjustmentMultiplicative(Finance::TermStructure<double, double> & sSigmaOISTS, 
                                                      Finance::TermStructure<double, double> & sSigmaCollatTS, 
                                                      const double dLambdaOIS, 
                                                      const double dLambdaCollat, 
                                                      const double dRhoCollatOIS,
                                                      const double dt,
                                                      const double dT1,
                                                      const double dT2,
                                                      const std::size_t iNIntervals) const;
	
		virtual double SwapQuantoAdjustmentMultiplicative(const Finance::TermStructure<double, double> & sSigmaOISTS, 
													   const Finance::TermStructure<double, double> & sSigmaCollatTS,
                                                       double dLambdaOIS, 
                                                       double dLambdaCollat, 
                                                       double dRhoCollatOIS,
													   const Finance::YieldCurve & sYieldCurveOIS,
													   const Finance::YieldCurve & sYieldCurveCollat,
                                                       double dt,
													   const std::vector <double> dS,
													   const std::vector <double> dT) const;
	
		};
}

#endif
