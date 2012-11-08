//
//  Test.h
//  Seminaire
//
//  Created by Alexandre HUMEAU on 06/11/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef Seminaire_Test_h
#define Seminaire_Test_h

#include <vector>

namespace Test {
    
    class TestVectorMother
    {
    public:
        std::vector<double> dVectorTest_;
        std::size_t i_;
    
        TestVectorMother(const std::vector<double> & dVectorTest);
        ~TestVectorMother();
        
        virtual void DoSomething();
    };
    
    class TestVectorDaughter : public TestVectorMother
    {
    public:
        TestVectorDaughter(const std::vector<double> & dVectorTest);
        ~TestVectorDaughter();
        
        virtual void DoSomething();
    };
}

#endif
