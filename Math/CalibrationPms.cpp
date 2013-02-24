//
//  CalibrationPms.cpp
//  Seminaire
//
//  Created by Alexandre HUMEAU on 23/02/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include "CalibrationPms.h"

namespace Calibration {
    
    CalibrationPms::CalibrationPms()
    {
        //  initialize value for Newton-Raphson Algorithm
        dLambda_ = 0.05;
        dSigma_ = 0.01;
    }
    
    CalibrationPms::~CalibrationPms()
    {}
    
    double CalibrationPms::GetLambda() const
    {
        return dLambda_;
    }
    
    double CalibrationPms::GetSigma() const
    {
        return dSigma_;
    }
    
    std::vector<double> CalibrationPms::LoadDataFromFile(const std::string &cFileName) const
    {
        std::ifstream sFile;
        sFile.open(cFileName.c_str());
        if (sFile.is_open())
        {
            std::vector<double> dResults;
            std::string cLine;
            while (std::getline(sFile, cLine)) 
            {
                std::istringstream iss(cLine);
                double dValue;
                while (iss >> dValue)
                {
                    dResults.push_back(dValue);
                }
            }
            return dResults;
        }
        else 
        {
            throw "Could not open file : " + cFileName;
        }
    }
    
}