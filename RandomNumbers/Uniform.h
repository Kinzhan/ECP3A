//
//  Uniform.h
//  MyLibrary
//
//  Created by Alexandre HUMEAU on 27/06/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef MyLibrary_Uniform_h
#define MyLibrary_Uniform_h

#include <vector>

namespace RandomNumbers {
    
    // The following class is used to generate uniform random variables 
    class Uniform
    {
    public:
        Uniform();
        Uniform(double dLeft, double dRight, std::size_t iNRealisations, int iAntitheticVariables = false);
        virtual ~Uniform();

        virtual void GenerateUniform();
        virtual std::vector<double> GetRealisations() const
        {
            return dRealisations_;
        };

    protected:
        double dLeft_;
        double dRight_;
        int iAntitheticVariables_;

        std::vector<double> dRealisations_;
        size_t iNRealisations_;

    };

}
#endif
