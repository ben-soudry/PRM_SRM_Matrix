#include "prm_srm_dot_prod.hpp"


int main(){
    //std::string fastaFile = "uniprot-proteome_homo_sapiens_reference_proteome_5_23_2016.fasta";

    std::string prmFile = "result/PRM_20150306_Rix_Ceri-heavy_PolyMac_10a.mzid_20150306_Rix_Ceri-heavy_PolyMac_10a.mgf";

    std::string resultsTsvFilename = "result/search_20150306_Rix_Ceri-heavy_PolyMac_10a.mzid_20150306_Rix_Ceri-heavy_PolyMac_10a.mzid.tsv"; 

    PrmSrmMatrix prmSrmMatrix;

    std::string path = "./result";
    std::vector<std::string> filenames;

    for (const auto & entry : std::experimental::filesystem::directory_iterator(path)){
        filenames.push_back(entry.path());
    }
    
     
    std::vector<std::string> spectra_filenames;
    std::vector<std::string> tsv_filenames;
    
    for(int i = 0; i < filenames.size(); i++){
         std::string last = filenames[i].substr(filenames[i].length() - 4);

 
         if(last == ".mgf"){
            spectra_filenames.push_back(filenames[i]); 
         }
         else if(last == ".tsv"){
            tsv_filenames.push_back(filenames[i]);   
         }
    } 
 
    //Run processSpectra on all the pairs of spectra/tsv files to build the matrix
        
    for(int i = 0; i < spectra_filenames.size(); i++){
        for(int j = 0; j < tsv_filenames.size(); j++){

                                    
            std::string tsv_basename = tsv_filenames[j].substr(tsv_filenames[j].find_last_of("/\\") + 1);

            std::string spectra_basename = spectra_filenames[i].substr(spectra_filenames[i].find_last_of("/\\") + 1);
                        
            spectra_basename = spectra_basename.substr(0, spectra_basename.length() - 4);
            
            //std::cout << "processingSpectra on: "  << tsv_basename << " and " << spectra_basename << std::endl;

            if(tsv_basename == "search_" + spectra_basename + ".mzid.tsv"){
                if(tsv_filenames[j] != "./result/search_20150306_Rix_Ceri-heavy_PolyMac_10a.mzid_20150306_Rix_Ceri-heavy_PolyMac_10a.MGF.mzid.tsv") {
                    continue;
                }

                std::cout << "=========================================================================" << std::endl; 
                std::cout << "processingSpectra on: "  << spectra_filenames[i] << std::endl;
                std::cout << " and TSV: " << tsv_filenames[j] << std::endl;
                
                auto startTime = std::chrono::system_clock::now();
                prmSrmMatrix.processSpectra(spectra_filenames[i], tsv_filenames[j]);
                    
                auto endTime = std::chrono::system_clock::now();
                std::chrono::duration<double> elapsed_seconds = endTime-startTime;
                std::cout << " Finished processingSpectra in " << elapsed_seconds.count() << " seconds." << std::endl;
             }
        }
    }
   }


PrmSrmMatrix::PrmSrmMatrix(){
    /*ranges = {-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0,
               0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,
               10.0,11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
               20.0,21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0};
    ranges_prm_plus_srm
           = {-18.0, -16.0, -14.0, -12.0, -10.0, -8.0, -6.0, -4.0, -2.0,
               0.0,  2.0,  4.0,  6.0,  8.0,  10.0,  12.0,  14.0,  16.0,  18.0,
               20.0,22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0,
               40.0,42.0, 44.0, 46.0, 48.0, 50.0, 52.0, 54.0, 58.0, 60.0};
      */ 
   
    ranges = {-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0,
               0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,
               10.0,11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0};
    ranges_prm_plus_srm
           = {-18.0, -16.0, -14.0, -12.0, -10.0, -8.0, -6.0, -4.0, -2.0,
               0.0,  2.0,  4.0,  6.0,  8.0,  10.0,  12.0,  14.0,  16.0,  18.0,
               20.0,22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0};
   
  
    for(int i = 0; i < numRanges; i++){
        prm_countMatrix[0][i] = 0;
        prm_countMatrix[1][i] = 0;
        srm_countMatrix[0][i] = 0;
        srm_countMatrix[1][i] = 0;
        prm_plus_srm_countMatrix[0][i] = 0;
        prm_plus_srm_countMatrix[1][i] = 0;
    } 
    totalCounts = 0;
}


void PrmSrmMatrix::processSpectra(std::string prmSrmFile, std::string resultsTsvFile)
{
    std::cout << "Getting experimental PRM-SRM" << std::endl;
    std::cout << "=================================================" << std::endl;
    
    //Get experimental PRM SRM from specta
    auto spectraPrmSrm = new prmSrm(prmSrmFile);
    
    modify_spectra(*spectraPrmSrm);

    std::unordered_map<int, experimentalPrmSrm> experimentalPrmSrmMap;

    for(int i = 0; i < spectraPrmSrm->n_spectra; i++){
        experimentalPrmSrm e;
        std::cout << "Length: " << spectraPrmSrm->spectra_length[i] << std::endl;
     
        auto e_prm_vec = std::vector<double>(spectraPrmSrm->prm[i], 
                                            spectraPrmSrm->prm[i] + spectraPrmSrm->spectra_length[i]);
         
        auto e_srm_vec = std::vector<double>(spectraPrmSrm->srm[i], 
                                            spectraPrmSrm->srm[i] + spectraPrmSrm->spectra_length[i]);
        
        //Make combined PRM+SRM vector
        std::vector<double> e_prm_plus_srm_vec;
        for(int i = 0; i < e_prm_vec.size(); i++){
            e_prm_plus_srm_vec.push_back(e_prm_vec[i]+e_srm_vec[i]);
        }
        
        e.prm = e_prm_vec;
        e.srm = e_srm_vec;
        e.prm_plus_srm = e_prm_plus_srm_vec;
        e.charges = spectraPrmSrm->charges[i];
        e.precursor_mass = spectraPrmSrm->precursor_mass[i];

        experimentalPrmSrmMap.insert(std::make_pair(spectraPrmSrm->spectra_index[i], e));
    }
    

    std::cout << "Getting MSGF+ spectra-peptide matching results" << std::endl;
    std::cout << "=================================================" << std::endl;
        
    auto tsvResults = new ParseResultsTsv(resultsTsvFile);
    
    std::cout << "Making PRM-SRM distribution matrix" << std::endl;
    std::cout << "=================================================" << std::endl;
    
    int matchCount = 0;
     
    std::unordered_map<int, bool> seen_specID;

    for(int i = 0; i < tsvResults->matches.size(); i++){
        //for(int j = 0; j < tsvResults->matches[i].proteins.size(); j++){

        std::string first = "FEIINAIY"; 
        std::string peptideFirst = tsvResults->matches[i].peptide.substr(0, 8);        


        if(tsvResults->matches[i].specID != 4404 ||  peptideFirst != first){
            continue;
        }

        if(tsvResults->matches[i].eValue >= 0.001 /*0.000000000001*/){
            continue;
        }       
        tsvResults->matches[i].peptide.erase(std::remove_if(tsvResults->matches[i].peptide.begin(), tsvResults->matches[i].peptide.end(), 
[]( auto const& c ) -> bool { return !std::isalnum(c); } ), tsvResults->matches[i].peptide.end()); 

        std::cout << "***Peptide: " << tsvResults->matches[i].peptide << " SpecID: "<< tsvResults->matches[i].specID << std::endl; 

        double* t_prm = nullptr; //Note, theoretical_vector allocates the memory
        double* t_srm = nullptr;
        double total_mass;
        theoretical_vector(0, tsvResults->matches[i].peptide, t_prm, t_srm, total_mass);

        auto t_prm_vec = std::vector<double>(t_prm, t_prm + tsvResults->matches[i].peptide.length()-1);
        auto t_srm_vec = std::vector<double>(t_srm, t_srm + tsvResults->matches[i].peptide.length()-1);
                 
        //std::cout << "Matching protein: " << tsvResults->matches[i].title << std::endl;
        //std::cout << " SpecID: " <<  tsvResults->matches[i].specID << std::endl;
        //std::cout << " E-Value: " << tsvResults->matches[i].eValue << std::endl;


        //Put peaks in a hash table for fast access
       
        std::unordered_map<int, bool> has_t_prm_peak;
        std::unordered_map<int, bool> has_t_srm_peak;
        
        std::cout << " PRM: "; 
        for(int k = 0; k < t_prm_vec.size(); k++){
            std::cout << t_prm_vec[k] << ", ";
            has_t_prm_peak.insert(std::make_pair((int) std::round(t_prm_vec[k]), true));
        }
        std::cout << std::endl;
        std::cout << " SRM: "; 
        for(int k = 0; k < t_srm_vec.size(); k++){
            std::cout << t_srm_vec[k] << ", ";
            has_t_srm_peak.insert(std::make_pair((int) std::round(t_srm_vec[k]), true));
        }
        std::cout << std::endl;  
        //Filling count matrix
        int specID =  tsvResults->matches[i].specID;
        auto emapFind = experimentalPrmSrmMap.find(specID);
        if(emapFind != experimentalPrmSrmMap.end()){
            matchCount += 1;
            //std::cout << "Filling count matrix!" << std::endl;
            experimentalPrmSrm e_prm_srm  = emapFind->second;

            
            double result = getDotProd(e_prm_srm, t_prm_vec, t_srm_vec);

    
            std::cout << "*************************** " << std::endl;
            std::cout << "*****Result! = " << result << std::endl;
            std::cout << "*************************** " << std::endl;
            std::cout << "Length " << e_prm_srm.prm.size() << std::endl;
            for(int i = 0; i < e_prm_srm.prm.size(); i++){
                totalCounts += 1;
                int prm_rangeId = getRangeID(e_prm_srm.prm[i]);
                 
                int prm_plus_srm_rangeId = getRangePrmSrmID(e_prm_srm.prm_plus_srm[i]);
 
                //std::cout << "prm_range_id: " << prm_rangeId << std::endl;
                auto has_prm_peak = has_t_prm_peak.find(i);
                if(has_prm_peak != has_t_prm_peak.end()){
                    prm_countMatrix[1][prm_rangeId] += 1;
                    prm_plus_srm_countMatrix[1][prm_rangeId] += 1;
                } else {
                    prm_countMatrix[0][prm_rangeId] += 1;
                    prm_plus_srm_countMatrix[0][prm_rangeId] += 1;
                }
                
                int srm_rangeId = getRangeID(e_prm_srm.srm[i]);
                //std::cout << "srm_range_id: " << srm_rangeId << std::endl;
                
                auto has_srm_peak = has_t_srm_peak.find(i);
                if(has_srm_peak != has_t_srm_peak.end()){
                    srm_countMatrix[1][srm_rangeId] += 1;
                    prm_plus_srm_countMatrix[1][srm_rangeId] += 1;
                } else {
                    srm_countMatrix[0][srm_rangeId] += 1;
                    prm_plus_srm_countMatrix[0][srm_rangeId] += 1;
                }
                
                
                //PRM+SRM combined
                /*if(has_prm_peak != has_t_prm_peak.end() || has_srm_peak != has_t_srm_peak.end()  ){
                    prm_plus_srm_countMatrix[1][prm_plus_srm_rangeId] += 1;
                } else {
                    prm_plus_srm_countMatrix[0][prm_plus_srm_rangeId] += 1;
                }*/
            }

        }                
        //std::cout << std::endl; 
    }
    delete spectraPrmSrm;
    delete tsvResults;

    
    std::cout << "# Matches Found: " << matchCount << std::endl;

    
     
    /*std::cout << "Results: " << std::endl;
    std::cout << "=========================================================================" << std::endl; 
    std::cout << "Total Matches considered: " << totalCounts << std::endl;
     
    std::cout << "PRM Count Matrix: " << std::endl;
    for(int i = 0; i < numRanges; i++){
        std::cout << prm_countMatrix[0][i] << ", ";
    }
    std::cout << std::endl << std::endl;
    for(int i = 0; i < numRanges; i++){
        std::cout << prm_countMatrix[1][i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "=================================================" << std::endl;
    
    std::cout << "SRM Count Matrix: " << std::endl;
    for(int i = 0; i < numRanges; i++){
        std::cout << srm_countMatrix[0][i] << ", ";
    }
    std::cout << std::endl << std::endl;
    for(int i = 0; i < numRanges; i++){
        std::cout << srm_countMatrix[1][i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "=================================================" << std::endl;
    
    std::cout << "PRM+SRM Count Matrix: " << std::endl;
    for(int i = 0; i < numRanges; i++){
        std::cout << prm_plus_srm_countMatrix[0][i] << ", ";
    }
    std::cout << std::endl << std::endl;
    for(int i = 0; i < numRanges; i++){
        std::cout << prm_plus_srm_countMatrix[1][i] << ", ";
    }
     
    std::cout << std::endl;
    std::cout << "=================================================" << std::endl;

    std::cout << "PRM prob Matrix: " << std::endl;
    for(int i = 0; i < numRanges; i++){
        prm_probMatrix[0][i] = (double) prm_countMatrix[0][i] / (double) totalCounts;
        std::cout << prm_probMatrix[0][i] << ", ";
    }
    std::cout << std::endl << std::endl;
    for(int i = 0; i < numRanges; i++){
        prm_probMatrix[1][i] = (double) prm_countMatrix[1][i] / (double) totalCounts;
        std::cout << prm_probMatrix[1][i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "=================================================" << std::endl;
 
    std::cout << "SRM prob Matrix: " << std::endl;
    for(int i = 0; i < numRanges; i++){
        srm_probMatrix[0][i] = (double) srm_countMatrix[0][i] / (double) totalCounts;
        std::cout << srm_probMatrix[0][i] << ", ";
    }
    std::cout << std::endl << std::endl;
    for(int i = 0; i < numRanges; i++){
        srm_probMatrix[1][i] = (double) srm_countMatrix[1][i] / (double) totalCounts;
        std::cout << srm_probMatrix[1][i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "=================================================" << std::endl;
     
    std::cout << "PRM+SRM prob Matrix: " << std::endl;
    for(int i = 0; i < numRanges; i++){
        prm_plus_srm_probMatrix[0][i] = (double) prm_plus_srm_countMatrix[0][i] / (double) 2*totalCounts;
        std::cout << prm_plus_srm_probMatrix[0][i] << ", ";
    }
    std::cout << std::endl << std::endl;
    for(int i = 0; i < numRanges; i++){
        prm_plus_srm_probMatrix[1][i] = (double) prm_plus_srm_countMatrix[1][i] / (double) 2*totalCounts;
        std::cout << prm_plus_srm_probMatrix[1][i] << ", ";
    }
    std::cout << std::endl;*/
    std::cout << "=================================================" << std::endl;
}


int PrmSrmMatrix::getRangeID(double peakMagnitude){
    for(int i = 0; i < ranges.size(); i++){
        if(peakMagnitude < ranges[i]){
            return i;
        }
    }
    //Larger than all ranges, return last index.
    return this->numRanges-1;

}

int PrmSrmMatrix::getRangePrmSrmID(double peakMagnitude){
    for(int i = 0; i < ranges_prm_plus_srm.size(); i++){
        if(peakMagnitude < ranges_prm_plus_srm[i]){
            return i;
        }
    }
    //Larger than all ranges, return last index.
    return this->numRanges-1;

}


double PrmSrmMatrix::getDotProd(experimentalPrmSrm e, std::vector<double> t_prm, std::vector<double> t_srm){

      std::map<int, bool> has_t_prm_peak;
      std::map<int, bool> has_t_srm_peak;
      
 
      //std::cout << " PRM: "; 
      for(int k = 0; k < t_prm.size(); k++){
        //std::cout << t_prm_vec[k] << ", ";
        has_t_prm_peak.insert(std::make_pair((int) std::round(t_prm[k]), true));
      }
      //std::cout << std::endl;
      //std::cout << " SRM: "; 
      for(int k = 0; k < t_srm.size(); k++){
        //std::cout << t_srm_vec[k] << ", ";
        has_t_srm_peak.insert(std::make_pair((int) std::round(t_srm[k]), true));
      }
      double result = 0.0;
      for(auto it = has_t_prm_peak.begin(); it != has_t_prm_peak.end(); ++it){


        if((*it).first < e.prm.size()){
            std::cout << "PRM peak at " <<  (*it).first << " Magnitude: " << e.prm[(*it).first] << std::endl;
            result += e.prm.at((*it).first);
        }
        else {
            std::cout << "PRM peak at " <<  (*it).first << " out of bounds " << std::endl;
        }
      }

      for(auto it = has_t_srm_peak.begin(); it != has_t_srm_peak.end(); ++it){
        if((*it).first < e.srm.size()){
            std::cout << "SRM peak at " <<  (*it).first << " Magnitude: " << e.srm[(*it).first] << std::endl;
            result += e.srm.at((*it).first);
        }
        else {
            std::cout << "SRM peak at " <<  (*it).first << " out of bounds " << std::endl;
        }
    
      }
      std::cout << "spectra print: " << std::endl;
      for(int i = 0; i < e.prm.size();  i++){
        std::cout << e.prm[i] << ", ";
      }
      for(int i = 0; i < e.srm.size();  i++){
        std::cout << e.srm[i] << ", ";
      }
      std::cout << std::endl;

      return result;
 

}


void PrmSrmMatrix::modify_spectra(prmSrm& ps) {
    for (int i = 0; i < ps.n_spectra; i++) {
        int len = ps.spectra_length[i];
       
        std::map <double, std::pair<int, bool> > bigprmsrm;

        for (int j = 0; j < len; j++) {
            if (ps.prm[i][j] < 1) ps.prm[i][j] = 0;
            //else ps.prm[i][j] = sqrt(ps.prm[i][j]);
            if (ps.srm[i][j] < 1) ps.srm[i][j] = 0;
            //else ps.srm[i][j] = sqrt(ps.srm[i][j]);
            
            if (bigprmsrm.size() < 20) bigprmsrm.insert(std::make_pair(ps.prm[i][j], std::make_pair(j, 0)));
            else if (bigprmsrm.begin()->first < ps.prm[i][j]) {
                bigprmsrm.erase(bigprmsrm.begin());
                bigprmsrm.insert(std::make_pair(ps.prm[i][j], std::make_pair(j, 0)));
            }

            if (bigprmsrm.size() < 20) bigprmsrm.insert(make_pair(ps.srm[i][j], std::make_pair(j, 1)));
            else if (bigprmsrm.begin()->first < ps.srm[i][j]) {
                bigprmsrm.erase(bigprmsrm.begin());
                bigprmsrm.insert(std::make_pair(ps.srm[i][j], std::make_pair(j, 1)));
            }
        }

        for (auto it = bigprmsrm.begin(); it != bigprmsrm.end(); ++it) {
            int idx = it->second.first;
            bool srmyes = it->second.second;
            if (srmyes) ps.srm[i][idx] /= 10000;
            else ps.prm[i][idx] /= 10000;
        }
    }
}




