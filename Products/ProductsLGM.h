//
//  ProductsLGM.h
//  Seminaire
//
//  Created by Alexandre HUMEAU on 31/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef Seminaire_ProductsLGM_h
#define Seminaire_ProductsLGM_h

#include <iostream>
#include "HullWhite.h"

namespace Products {
    class ProductsLGM : public Processes::LinearGaussianMarkov
    {
    public:
        
        ProductsLGM(const Processes::LinearGaussianMarkov & sLGMProcess, const double dEpsilonMaturity = 0.001);
        
        virtual std::vector<double> RiskNeutralDiscountFactor(const std::size_t iPath, const Finance::SimulationData & sSimulationData) const;
        
        virtual std::vector<double> Caplet(const double dStart, const double dEnd, const double dStrike, const Finance::SimulationData & sSimulationData) const;
    private:
        //  about 1 day
        double dEpsilonMaturity_;
    };
}

#endif
