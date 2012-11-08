//
//  HJM.h
//  Seminaire
//
//  Created by Alexandre HUMEAU on 24/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef Seminaire_HJM_h
#define Seminaire_HJM_h

//////////////////////////////////////////////////////////////////////////////////
//  
//  We are considering here a Hull-White diffusion of the short rate.
//  But in fact, we consider a HJM framework to diffuse the 
//  instantaneous forward rate with a volatility, which makes the
//  IFR a Markov process
//
//  For now, we choose to simulate only forward rates, which are given
//  in the initial yield curve. We can in a second step, put pillars
//  in the yield curve and choose to simulate these pillars, which are
//  found for date 0, with a proper spline cubic interpolation
//
//  The IFR follows the following SDE
//  df(t,T) = \sigma(t,T) \Sigma(t,T) dt + \sigma(t,T) dW_t
//
//  If we are in a Hull-White case, \sigma(t,T) = \sigma_t exp(-\lambda (T-t))
//  \Simga(t,T) = \int_t^T \sigma(s,T) ds
//
/////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include "YieldCurve.h"
#include "SimulationData.h"

namespace Processes {
    
    enum HJMType
    {
        HULL_WHITE
    };
    
    enum SimulationProbability
    {
        T_FORWARD_NEUTRAL,
        RISK_NEUTRAL
    };
    
    class HeathJarrowMorton
    {
    protected:
        //  Initial YieldCurve
        Finance::YieldCurve sInitialYieldCurve_;
        
    public:
        HeathJarrowMorton();
        HeathJarrowMorton(const Finance::YieldCurve & sInitialYieldCurve);
        ~HeathJarrowMorton();
        
        virtual void Simulate(const std::size_t iNRealisations,
                              const std::vector<double> & dSimulationTenors,
                              Finance::SimulationData & sSimulationData,
                              bool bIsStepByStepMC) const;
        virtual double A(double t) const = 0;
        virtual double B(double t) const = 0;
        virtual double BondPrice(const double dt, const double dMaturity, const double dX, const SimulationProbability eProbability) const;
        virtual double Libor(const double dt, const double dStart, const double dEnd, const double dX, const SimulationProbability eProbability) const;
    };
    
}
#endif
