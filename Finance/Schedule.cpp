//
//  Schedule.cpp
//  Seminaire
//
//  Created by Alexandre HUMEAU on 23/01/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "Schedule.h"

namespace Finance {
    
    Schedule::Schedule(const Utilities::Date::MyDate & sStart, const Utilities::Date::MyDate & sEnd, const YieldCurve & sYieldCurve, const MyBasis eBasis, const MyFrequency eFrequency) : eFrequency_(eFrequency)
    {
        Utilities::Date::MyDate sCurrentStart, sCurrentEnd = sEnd;
        std::pair<std::size_t, Utilities::Date::TimeUnits> NumberAndUnitToAdd = Frequency::ParseFrequency(eFrequency_);
        sCurrentStart = sCurrentEnd;
        sCurrentEnd.Add( - NumberAndUnitToAdd.first, NumberAndUnitToAdd.second);
        
        while (sCurrentStart > sStart)
        {
            EventOfSchedule sEvent(sCurrentStart, sCurrentEnd, sYieldCurve, eBasis);
            sSchedule_.push_back(sEvent);
            
            //  Update current
            sCurrentEnd.Add( - NumberAndUnitToAdd.first, NumberAndUnitToAdd.second);
        }
    }
    
    Schedule::~Schedule()
    {
        sSchedule_.~vector();
    }
    
    std::vector<EventOfSchedule> Schedule::GetSchedule() const
    {
        return sSchedule_;
    }
}