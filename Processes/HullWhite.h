//
//  HullWhite.h
//  Seminaire
//
//  Created by Alexandre HUMEAU on 19/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef Seminaire_HullWhite_h
#define Seminaire_HullWhite_h

#include <iostream> 
#include "HJM.h"

namespace Processes {

    class LinearGaussianMarkov : public HeathJarrowMorton
    {
    protected:
        //  Parameters of the model
        double dLambda_;
        Finance::TermStructure<double, double> dSigma_;
        
    public:
		
        LinearGaussianMarkov();
        LinearGaussianMarkov(const Finance::YieldCurve & sDiscountCurve, double dLambda, const Finance::TermStructure<double, double> & dSigma);
        LinearGaussianMarkov(const Finance::YieldCurve & sDiscountCurve, const Finance::YieldCurve & sForwardCurve, double dLambda, const Finance::TermStructure<double, double> & dSigma);
 
        virtual ~LinearGaussianMarkov();

        virtual double GetLambda() const
        {
            return dLambda_;
        }
        
        virtual Finance::TermStructure<double, double> GetSigma() const
        {
            return dSigma_;
        }
        
        virtual void SetSigma(const Finance::TermStructure<double, double> & sSigmaTS)
        {
            dSigma_ = sSigmaTS;
        }
        
        virtual double BondPrice(double dt, double dT, double dX, const CurveName & eCurveName) const;
        virtual double Libor(double dt, double dStart, double dEnd, double dX, const CurveName & eCurveName, double dQA = 1.0) const;
        virtual void Simulate(std::size_t iNRealisations,
                              const std::vector<double> & dSimulationTenors,
                              Finance::SimulationData & sSimulationData,
                              bool bIsStepByStepMC) const;
        virtual double BracketChangeOfProbability(double dt, double dT) const;
        virtual void ChangeOfProbability(double dT, const Finance::SimulationData & sSimulationDataRiskNeutral,
                                         Finance::SimulationData & sSimulationDataTForward) const;
        
        virtual double DeterministPart(double dt, double dT) const;
        virtual double A(double t) const;
        virtual double B(double t) const;
    };
}

#endif
