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

namespace Finance
{
    Annuity::~Annuity()
    {
        sStart_.~MyDate();
        sEnd_.~MyDate();
        sYieldCurve_.~YieldCurve();
    }
    
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
        //  Get the schedule
        Schedule sSchedule(sStart_, sEnd_, sYieldCurve_, eBasis_, eFrequency_);
        
        //  Definition of annuity : 
        //  Level = \sum_{i} (\tau(T_i, T_{i+1}) * DF(t,T_{i+1}))
        //  where \tau(T_i, T_{i+1}) is the coverage (year fraction) between T_i, T_{i+1} and DF is the discount factor function.
        
        double dLevel = 0;
        for (std::size_t iDate = 0 ; iDate < sSchedule.GetSchedule().size() ; ++iDate)
        {
            dLevel += sSchedule.GetSchedule()[iDate].GetCoverage() * sSchedule.GetSchedule()[iDate].GetPayingDateDF();
        }
        
        return dLevel;
    }
}