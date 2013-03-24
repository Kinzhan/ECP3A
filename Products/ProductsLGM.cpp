//
//  ProductsLGM.cpp
//  Seminaire
//
//  Created by Alexandre HUMEAU on 31/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "ProductsLGM.h"

namespace Products {
    
    ProductsLGM::~ProductsLGM()
    {}
    
    ProductsLGM::ProductsLGM(const Processes::LinearGaussianMarkov & sLGMProcess, double dEpsilonMaturity) : dEpsilonMaturity_(dEpsilonMaturity), LinearGaussianMarkov(sLGMProcess)
    {}
    
    std::vector<double> ProductsLGM::Caplet(double dStart, double dEnd, double dPay, double dStrike, const Finance::SimulationData & sSimulationData, const Processes::CurveName & eCurveName, double dQA) const
    {
        //  Price of a caplet starting a dStart, ending at dEnd and paying at dPay, with Strike dStrike and with MC Simulation factors at dStart
        std::size_t iWhere = Utilities::FindInVector(sSimulationData.GetDateList(), static_cast<long>(dStart * 365));
        
        if (iWhere != sSimulationData.GetDateList().size())
        {
            std::vector<std::vector<double> > dMatrixEndFactor = sSimulationData.GetData().second[iWhere];
        
            std::size_t iNPaths = dMatrixEndFactor.size();
            std::vector<double> dResults;
        
            double dCoverage = (dEnd - dStart);
            
            for (std::size_t iPath = 0 ; iPath < iNPaths ; ++iPath)
            {
                //  Only one factor which is simulated for now
                //  Fixing of the libor at start date of the period
                //  Alexandre 4/12/2012 add coverage because cash-flow of cash-flow is cvg * max (Libor - K, 0)
                double dFactor = dMatrixEndFactor[iPath][0];
                double dLibor = Libor(dStart, dStart, dEnd, dFactor/*, Processes::T_FORWARD_NEUTRAL*/, eCurveName, dQA);
                dResults.push_back( dCoverage * std::max(dLibor - dStrike, 0.0));
            }
        
            return dResults;
        }
        else
        {
            std::cout<< "Start date not found in simulation" << std::endl;
            return std::vector<double>(0,0.0);
        }
    }
    
    std::vector<double> ProductsLGM::RiskNeutralDiscountFactor(const std::size_t iPath, const Finance::SimulationData &sSimulationData) const
    {
        // not implemented
        return std::vector<double>(0,0.0);
    }
    
}