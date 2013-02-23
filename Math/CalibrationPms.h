//
//  CalibrationPms.h
//  Seminaire
//
//  Created by Alexandre HUMEAU on 23/02/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#ifndef Seminaire_CalibrationPms_h
#define Seminaire_CalibrationPms_h

#include <vector>

namespace Calibration {
    class CalibrationPms
    {
    protected:
        double dLambda_, dSigma_;
    public:
        CalibrationPms();
        virtual ~CalibrationPms();
        
        //  Getters
        virtual double GetLambda() const;
        virtual double GetSigma() const;
        
        virtual std::vector<double> LoadDataFromFile(const std::string & cFileName) const;
    };
}

#endif
