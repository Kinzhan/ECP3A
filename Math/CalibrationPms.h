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
    
    //  Maximum Likelihood estimation of the paramaters knowing the historical fixing of the diffusion
    //  We assume a Ornstein Ulhenbeck process for the Instantaneous Forward Rate
    //  dX_t = \lambda (\mu - X_t) dt + \sigma dW_t
    
    struct NewtonPms
    {
        NewtonPms(double dTolerance, std::size_t iNIterMax);
        
        double dTolerance_; // tolerance for lambda in the Newton
        std::size_t iNIterMax_; // number of iterations in the Newton
    };
    
    struct NewtonFunctionMu0
    {
        double dDeltaT_;
        std::vector<double> dData_;
        double dEpsilonDerivative_;
        NewtonFunctionMu0(double dDeltaT, const std::vector<double> & dData, double dEpsilonDerivative = 0.00001);
        
        double func(double dLambda);
        double dfdlambda(double dLambda);
    };
    
    struct NewtonFunction
    {
        double dDeltaT_;
        std::vector<double> dData_;
        double dEpsilonDerivative_;
        NewtonFunction(double dDeltaT, const std::vector<double> & dData, double dEpsilonDerivative = 0.00001);
        
        double VolatilitySquare(double dLambda);
        double Mu(double dLambda);
        double func(double dLambda);
        double dfdlambda(double dlambda);
    };
    
    class CalibrationPms
    {
    protected:
        double dLambda_, dSigma_, dMu_;
        double dEpsilonDerivative_;
    public:
        CalibrationPms(double dLambda, double dEpsilonDerivative);
        virtual ~CalibrationPms();
        
        //  Getters
        virtual double GetLambda() const;
        virtual double GetSigma() const;
        virtual double GetMu() const;
        
        virtual std::vector<double> LoadDataFromFile(const std::string & cFileName, bool bIsLibor = false, double dDeltaT = 1./252.) const;
        
        virtual void ComputeSigmaMu0(const std::vector<double> & dData, double dDeltaT);
        virtual void ComputeSigma(const std::vector<double> & dData, double dDeltaT);
        virtual void ComputeMu(const std::vector<double> & dData, double dDeltaT);
        
        virtual void NewtonRaphsonAlgorithmMu0(const std::vector<double> & dData, double dDeltaT, const NewtonPms & sNewtonPms);
        virtual void NewtonRaphsonAlgorithm(const std::vector<double> & dData, double dDeltaT, const NewtonPms & sNewtonPms);
    };
}

#endif
