//
//  Annuity.cpp
//  Seminaire
//
//  Created by Alexandre HUMEAU on 23/01/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "Annuity.h"
#include "Schedule.h"
#include "DiscountFactor.h"
#include "Coverage.h"

namespace Finance
{
    Annuity::~Annuity()
    {}
    
    Annuity::Annuity(const Utilities::Date::MyDate     sStart,
                     const Utilities::Date::MyDate     sEnd,
                     const MyBasis                     eBasis,
                     const MyFrequency                 eFrequency,
                     const YieldCurve                  sYieldCurve): 
    
    sStart_(sStart), 
    sEnd_(sEnd), 
    eBasis_(eBasis), 
    eFrequency_(eFrequency),
    sYieldCurve_(sYieldCurve)
    
    {}
    
    double Annuity::ComputeAnnuity() const
    {
        DF sDF(sYieldCurve_);
        //  Get the schedule
        Schedule sSchedule(sStart_, sEnd_, sYieldCurve_, eBasis_, eFrequency_);
        
        //  Definition of annuity : 
        //  Level = \sum_{i} (\tau(T_i, T_{i+1}) * DF(t,T_{i+1}))
        //  where \tau(T_i, T_{i+1}) is the coverage (year fraction) between T_i, T_{i+1} and DF is the discount factor function.
        
        double dLevel = 0;
        std::vector<EventOfSchedule> sVectEventOfSchedule = sSchedule.GetSchedule();
        for (std::size_t iDate = 0 ; iDate < sVectEventOfSchedule.size() ; ++iDate)
        {
            dLevel += sVectEventOfSchedule[iDate].GetCoverage() * sVectEventOfSchedule[iDate].GetPayingDateDF();
        }
        //  add last date
        Coverage sCoverage(eBasis_, sVectEventOfSchedule.back().GetEndDate(), sEnd_);
        dLevel += sCoverage.ComputeCoverage() * sDF.DiscountFactor(sEnd_);
        
        return dLevel;
    }
}