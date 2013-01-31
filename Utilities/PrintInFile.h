//
//  PrintInFile.h
//  Seminaire
//
//  Created by Alexandre HUMEAU on 22/11/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef Seminaire_PrintInFile_h
#define Seminaire_PrintInFile_h

#include <vector>

namespace Utilities {
    
    class PrintInFile
    {
    public:
        PrintInFile();
        PrintInFile(const std::string & cFileName, const bool bAppend, const std::size_t iPrecision);
        virtual ~PrintInFile();
        
        virtual void PrintDataInFile(const std::vector<double> & dData);
        virtual void PrintDataInFile(const std::vector<std::pair<double, std::size_t> > & dData);
    protected:
        std::string cFileName_;
        bool bAppend_;
        std::size_t iPrecision_; // number of digits after the point for printing
        std::string cPrecision_; // Conversion of iPrecision_ into a std::string value
    };
}

#endif
