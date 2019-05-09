#include "parse_results_tsv.hpp"



ParseResultsTsv::ParseResultsTsv(std::string resultsTsvFilename){
    filename = resultsTsvFilename;
    std::cout << "  Parsing TSV..." << std::endl;
    std::ifstream tsvFile;

    tsvFile.open(filename);
    std::string colNames;
    std::getline( tsvFile, colNames);
    
    while(!tsvFile.eof()){
        Match match;

        std::string line;
        std::getline( tsvFile, line);
        if(line.length() == 0){
            continue;
        }

        std::stringstream buffer(line);
        std::string temp;
        
        getline( buffer, temp, '\t');
        match.specFile = temp;      
  
        getline( buffer, temp, '\t');
        //std::cout << temp.substr(6,100) << " | ";
        match.specID = std::stoi(temp.substr(6, 100));       

        getline( buffer, temp, '\t');
        //match.scanNum = std::stoi(temp);       

        getline( buffer, temp, '\t');
        match.title = temp;       

        getline( buffer, temp, '\t');
        match.fragMethod = temp;       
 
        getline( buffer, temp, '\t');
        match.fragMethod = temp;       

        getline( buffer, temp, '\t');
        match.precursor = temp;       

        getline( buffer, temp, '\t');
        match.precursorErrorDa = temp;       

        getline( buffer, temp, '\t');
        match.charge = temp;       

        getline( buffer, temp, '\t');
        match.peptide = temp;       

        getline( buffer, temp, '\t');
        //Get each matching protein
        std::stringstream proteinBuffer(temp);
        std::string p;
        while(proteinBuffer){
            getline(proteinBuffer, p, ';');
            match.proteins.push_back(p);
        }       
        
        getline( buffer, temp, '\t');
        match.deNovoScore = std::stod(temp);       
 
        getline( buffer, temp, '\t');
        match.MSGFScore = std::stod(temp);       
 
        getline( buffer, temp, '\t');
        match.specEValue = std::stod(temp);       
 
        getline( buffer, temp, '\t');
        match.eValue = std::stod(temp);       

        /*for(int i = 0; i < match.proteins.size(); i++){ 
            //std::cout << match.proteins[i] << ", ";
            
        }*/
        //std::cout << std::endl;
        
        matches.push_back(match);
    }
    


    tsvFile.close();
}
