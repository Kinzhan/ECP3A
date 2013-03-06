//
//  HJM.cpp
//  Seminaire
//
//  Created by Alexandre HUMEAU on 24/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "HJM.h"

namespace Processes {
    HeathJarrowMorton::HeathJarrowMorton()
    {}
    
    HeathJarrowMorton::HeathJarrowMorton(const Finance::YieldCurve & sDiscountCurve) : 
    sDiscountCurve_(sDiscountCurve)
    {}
    
    HeathJarrowMorton::HeathJarrowMorton(const Finance::YieldCurve & sDiscountCurve, const Finance::YieldCurve & sForwardCurve) : 
    sDiscountCurve_(sDiscountCurve)
    {
        sForwardCurve_ = sForwardCurve;
    }
    
    HeathJarrowMorton::~HeathJarrowMorton()
    {}
    
    /*void HeathJarrowMorton::Simulate(const std::size_t iNRealisations, 
                                     const std::vector<double> & dSimulationTenors ,//In years
                                     Finance::SimulationData & sSimulatedData,
                                     bool bIsStepbyStepMC) const
    {
        Utilities::require(!dSimulationTenors.empty(), "Simulation Times is empty");
        Utilities::require(iNRealisations > 0, "Number of paths has to be positive");
        
        //  In this function, we will simulate the factor \int_{0}^{t} a(s) dW^Q_s for t in the simulation tenors vector
        for (std::size_t iSimulationTenor = 0 ; iSimulationTenor < dSimulationTenors.size() ; ++iSimulationTenor)
        {
            //  We simulate the wanted factor
        }
    }*/
    
    /*double HeathJarrowMorton::BondPrice(const double dt, const double dMaturity, const double dX) const
    {
        return 0;
    }*/
    
    double HeathJarrowMorton::Libor(const double dt, const double dStart, const double dEnd, const double dX, const CurveName & eCurveName) const
    {
        //  Must change coverage to take into account real basis
        return 1.0 / (dEnd - dStart) * (BondPrice(dt, dStart, dX, eCurveName) / BondPrice(dt, dEnd, dX, eCurveName) - 1.0);
    }
}