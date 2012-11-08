//
//  Test.cpp
//  Seminaire
//
//  Created by Alexandre HUMEAU on 06/11/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "Test.h"

namespace Test {
    TestVectorMother::TestVectorMother(const std::vector<double> & dVector) : dVectorTest_(dVector), i_(0)
    {}
    
    TestVectorMother::~TestVectorMother()
    {
        //dVectorTest_.clear();
        //dVectorTest_.~vector();
    }
    
    void TestVectorMother::DoSomething()
    {
        i_++;
    }
    
    TestVectorDaughter::TestVectorDaughter(const std::vector<double > & dVector) : TestVectorMother(dVector)
    {}
    
    TestVectorDaughter::~TestVectorDaughter()
    {
        //TestVectorMother::~TestVectorMother();
        //dVectorTest_.clear();
        //dVectorTest_.~vector();
    }
    
    void TestVectorDaughter::DoSomething()
    {
        i_ += 2;
    }
    
    
}