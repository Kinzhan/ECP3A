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
    
    struct NewtonPms
    {
        NewtonPms(double dTolerance, std::size_t iNIterMax);
        
        double dTolerance_; // tolerance for lambda in the Newton
        std::size_t iNIterMax_; // number of iterations in the Newton
    };
    
    struct NewtonFunction
    {
        double dDeltaT_;
        std::vector<double> dData_;
        double dEpsilonDerivative_;
        NewtonFunction(double dDeltaT, const std::vector<double> & dData, double dEpsilonDerivative = 0.00001);
        
        double func(double dLambda);
        double dfdlambda(double dLambda);
    };
    
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
        
        virtual void ComputeSigma(const std::vector<double> & dData, double dDeltaT);
        
        virtual void NewtonRaphsonAlgorithm(const std::vector<double> & dData, double dDeltaT, const NewtonPms & sNewtonPms);
    };
}

#endif
