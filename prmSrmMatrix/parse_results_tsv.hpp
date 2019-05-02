#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>       

class ParseResultsTsv {

    public: 
    struct Match {
        std::string specFile;
        int specID;
        int scanNum;
        std::string title;
        std::string fragMethod;
        std::string precursor; 
        std::string precursorErrorDa;
        std::string charge;
        std::string peptide;
        std::vector<std::string> proteins;
        double deNovoScore;
        double MSGFScore;
        double specEValue;
        double eValue;
    };

    ParseResultsTsv(std::string resultsTsvFilename);
            


    std::vector<Match> matches;
    std::string filename;


};
