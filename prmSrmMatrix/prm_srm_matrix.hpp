#include "ppseq.hpp"
#include "io.hpp"
#include "theoretical_util.hpp"
#include "parse_results_tsv.hpp"


#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <chrono>
#include <array>
#include <vector>
#include <list>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <set>
#include <utility>
#include <numeric>
#include <algorithm>
#include <random>
#include <boost/functional/hash.hpp>
#include <unordered_map>


class PrmSrmMatrix {
    public:
        static constexpr int numRanges = 40;
       
        std::array<double, numRanges-1> ranges;
        PrmSrmMatrix(std::string fastaFile, std::string prmSrmFile, std::string resultsTsvFile);
        
        int getRangeID(double peakMagnitude);

        const int minPeptideLen = 10;
        const int maxPeptideLen = 50;
        
        struct theoreticalPrmSrm {
            std::vector<double> prm;
            std::vector<double> srm;
        };

         
        struct experimentalPrmSrm {
            std::vector<double> prm;
            std::vector<double> srm;
            int charges;
            double precursor_mass;
        };

        
        typedef std::array<std::array<int, numRanges>,2> CountMatrix;
        typedef std::array<std::array<double, numRanges>,2> ProbMatrix;

 
        CountMatrix prm_countMatrix;
        CountMatrix srm_countMatrix;
        int totalCounts;      
    
        ProbMatrix prm_probMatrix;
        ProbMatrix srm_probMatrix;

    private:        
        prmSrm* spectraPrmSrm;
        ppseq* peptideSeqDB;
        ParseResultsTsv* tsvResults; 

        int num_proteins;
        std::unordered_map<std::string, theoreticalPrmSrm> theoreticalPrmSrmMap;
        
};