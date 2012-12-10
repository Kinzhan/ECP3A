//
//  TermStructure.h
//  Seminaire
//
//  Created by Alexandre HUMEAU on 19/10/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef Seminaire_TermStructure_h
#define Seminaire_TermStructure_h

#include <vector>
#include "Require.h"

//  This file creates a termstructure template

namespace Finance {
    
    template<class T, class U>
    class TermStructure
    {
    public:
        TermStructure()
        {
            TVariables_.resize(1);
            TVariables_[0] = 0;
            
            UValues_.resize(1);
            UValues_[0] = 0;
        }
        
        TermStructure(const std::vector<T> & TVariables, const std::vector<U> & UValues) : TVariables_(TVariables), UValues_(UValues)
        {
            Utilities::require(TVariables.size() == UValues.size());
        }
        
        virtual ~TermStructure()
        {}
        
        virtual std::vector<T> GetVariables() const
        {
            return TVariables_;
        }
        
        virtual std::vector<U> GetValues() const
        {
            return UValues_;
        }
        
        virtual bool IsTermStructure() const
        {
            return (TVariables_.size() != 1) && (UValues_.size() != 1);
        }
        
        virtual U Interpolate(const T variable) const
        {
            //  Flat extrapolation on the left
            if (variable < TVariables_[0])
            {
                return UValues_[0];
            }
            for (std::size_t i = 0 ; i < TVariables_.size() - 1 ; ++i)
            {
                if (TVariables_[i] <= variable < TVariables_[i + 1] )
                {
                    return UValues_[i];
                }
            }
            
            //  Flat extrapolation on the right
            if (variable >= TVariables_.back())
            {
                return UValues_.back();
            }
            return 0.0;
        }
        
    private:
        std::vector<T> TVariables_;
        std::vector<U> UValues_;
    };
    
}

#endif
