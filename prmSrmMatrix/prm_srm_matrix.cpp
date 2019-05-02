#include "prm_srm_matrix.hpp"


int main(){

    std::string fastaFile = "uniprot-proteome_homo_sapiens_reference_proteome_5_23_2016.fasta";

    std::string prmFile = "result/PRM_20150306_Rix_Ceri-heavy_PolyMac_10a.mzid_20150306_Rix_Ceri-heavy_PolyMac_10a.mgf";

    std::string resultsTsvFilename = "result/search_20150306_Rix_Ceri-heavy_PolyMac_10a.mzid_20150306_Rix_Ceri-heavy_PolyMac_10a.mzid.tsv"; 

    PrmSrmMatrix prmSrmMatrix(fastaFile, prmFile, resultsTsvFilename);
}


PrmSrmMatrix::PrmSrmMatrix(std::string fastaFile, std::string prmSrmFile, std::string resultsTsvFile)
{
    ranges = {-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0,
               0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,
               10.0,11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
               20.0,21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0};
    
    for(int i = 0; i < numRanges; i++){
        prm_countMatrix[0][i] = 0;
        prm_countMatrix[1][i] = 0;
        srm_countMatrix[0][i] = 0;
        srm_countMatrix[1][i] = 0;
    } 
    totalCounts = 0;

    std::cout << "Getting experimental PRM-SRM" << std::endl;
    std::cout << "=================================================" << std::endl;
    
    //Get experimental PRM SRM from specta
    spectraPrmSrm = new prmSrm(prmSrmFile);
    
    std::unordered_map<int, experimentalPrmSrm> experimentalPrmSrmMap;

    for(int i = 0; i < spectraPrmSrm->n_spectra; i++){
        experimentalPrmSrm e;
        
        auto e_prm_vec = std::vector<double>(spectraPrmSrm->prm[i], 
                                            spectraPrmSrm->prm[i] + spectraPrmSrm->spectra_length[i]);
         
        auto e_srm_vec = std::vector<double>(spectraPrmSrm->srm[i], 
                                            spectraPrmSrm->srm[i] + spectraPrmSrm->spectra_length[i]);
        e.prm = e_prm_vec;
        e.srm = e_srm_vec;
        e.charges = spectraPrmSrm->charges[i];
        e.precursor_mass = spectraPrmSrm->precursor_mass[i];

        experimentalPrmSrmMap.insert(std::make_pair(spectraPrmSrm->spectra_index[i], e));
    }
    
    std::cout << "Getting peptides from protein database" << std::endl;
    std::cout << "=================================================" << std::endl;
    //Get peptide sequences from protein database
    read_ppDatabase(fastaFile, peptideSeqDB, num_proteins);
    

    /*for(int i = 0; i < num_proteins; i++){
        std::cout << "peptide: " <<   peptideSeqDB[i].seq << std::endl; 
    }*/

    std::cout << "Converting peptides to binary vectors" << std::endl;
    std::cout << "=================================================" << std::endl;
    
    //Convert peptides to binary vectors
    for(int i = 0; i < num_proteins; i++){
        //std::cout << peptideSeqDB[i].seq << std::endl; 
        if(peptideSeqDB[i].seq.length() >= minPeptideLen 
        && peptideSeqDB[i].seq.length() <= maxPeptideLen){
            //std::cout << "peptide: " << peptideSeqDB[i].name << std::endl; 
            double* t_prm = nullptr;
            double* t_srm = nullptr;
            double t_mass = 0;
            theoretical_vector(0, peptideSeqDB[i].seq, t_prm, t_srm, t_mass);
            auto t_prm_vec = std::vector<double>(t_prm, t_prm + peptideSeqDB[i].seq.length()-1);
            auto t_srm_vec = std::vector<double>(t_srm, t_srm + peptideSeqDB[i].seq.length()-1);
            
            theoreticalPrmSrm t;
            t.prm = t_prm_vec;
            t.srm = t_srm_vec; 
            
            //clean up protein name
            std::stringstream nameBuffer(peptideSeqDB[i].name);
            std::string proteinName;
            getline(nameBuffer, proteinName, ' ');

            std::cout << "peptide: " << proteinName << std::endl; 
            theoreticalPrmSrmMap.insert(std::make_pair(proteinName, t));
        }
    } 
    std::cout << "Getting MSGF+ spectra-peptide matching results" << std::endl;
    std::cout << "=================================================" << std::endl;
        
    tsvResults = new ParseResultsTsv(resultsTsvFile);

    
    std::cout << "Making PRM-SRM distribution matrix" << std::endl;
    std::cout << "=================================================" << std::endl;
    
    int matchCount = 0;
    for(int i = 0; i < tsvResults->matches.size(); i++){
        for(int j = 0; j < tsvResults->matches[i].proteins.size(); j++){
            if(tsvResults->matches[i].eValue >= 0.001){
                continue;
            }       
    
            std::string currProtein = tsvResults->matches[i].proteins[j]; 
            
            //clean up protein name
            std::stringstream nameBuffer(currProtein);
            std::string proteinName;
            getline(nameBuffer, proteinName, '(');

            
            std::cout << "Searching for protein: " << proteinName << std::endl;
            auto tmapFind = theoreticalPrmSrmMap.find(proteinName);

            if(tmapFind != theoreticalPrmSrmMap.end()){
                theoreticalPrmSrm t_prm_srm  = tmapFind->second;
                

                std::cout << "Found Match in Peptide database: " << proteinName;
                std::cout << " SpecID: " <<  tsvResults->matches[i].specID << std::endl;
                std::cout << "  E-Value: " << tsvResults->matches[i].eValue << std::endl;
                matchCount += 1;
                std::unordered_map<int, bool> has_t_prm_peak;
                std::unordered_map<int, bool> has_t_srm_peak;
                
                std::cout << " PRM: "; 
                for(int k = 0; k < t_prm_srm.prm.size(); k++){
                    std::cout << t_prm_srm.prm[k] << ", ";
                    has_t_prm_peak.insert(std::make_pair((int) std::round(t_prm_srm.prm[k]), true));
                }
                std::cout << std::endl;
                std::cout << " SRM: "; 
                for(int k = 0; k < t_prm_srm.srm.size(); k++){
                    std::cout << t_prm_srm.srm[k] << ", ";
                    has_t_srm_peak.insert(std::make_pair((int) std::round(t_prm_srm.srm[k]), true));
                }
                
                //Filling count matrix
                int specID =  tsvResults->matches[i].specID;
                auto emapFind = experimentalPrmSrmMap.find(specID);
                if(emapFind != experimentalPrmSrmMap.end()){
                    std::cout << "Filling count matrix!" << std::endl;
                    experimentalPrmSrm e_prm_srm  = emapFind->second;
                    for(int i = 0; i < e_prm_srm.prm.size(); i++){
                        totalCounts += 1;
                        int prm_rangeId = getRangeID(e_prm_srm.prm[i]);
                           
                        //std::cout << "prm_range_id: " << prm_rangeId << std::endl;
                        auto has_prm_peak = has_t_prm_peak.find(i);
                        if(has_prm_peak != has_t_prm_peak.end()){
                            prm_countMatrix[1][prm_rangeId] += 1;
                        } else {
                            prm_countMatrix[0][prm_rangeId] += 1;
                        }
                        
                        int srm_rangeId = getRangeID(e_prm_srm.srm[i]);
                        //std::cout << "srm_range_id: " << srm_rangeId << std::endl;
                        
                        auto has_srm_peak = has_t_srm_peak.find(i);
                        if(has_srm_peak != has_t_srm_peak.end()){
                            srm_countMatrix[1][srm_rangeId] += 1;
                        } else {
                            srm_countMatrix[0][srm_rangeId] += 1;
                        }
                    }

                }                
                
                std::cout << std::endl; 
            }
        }
    }
    std::cout << "# Matches Found: " << matchCount << std::endl;
    
    
    std::cout << "PRM prob Matrix: " << std::endl;
    for(int i = 0; i < numRanges; i++){
        prm_probMatrix[0][i] = (double) prm_countMatrix[0][i] / (double) totalCounts;
        std::cout << prm_probMatrix[0][i] << ", ";
    }
    std::cout << std::endl;
    for(int i = 0; i < numRanges; i++){
        prm_probMatrix[1][i] = (double) prm_countMatrix[1][i] / (double) totalCounts;
        std::cout << prm_probMatrix[1][i] << ", ";
    }
    std::cout << std::endl;
    
    std::cout << "SRM count Matrix: " << std::endl;
    for(int i = 0; i < numRanges; i++){
        srm_probMatrix[0][i] = (double) srm_countMatrix[0][i] / (double) totalCounts;
        std::cout << srm_probMatrix[0][i] << ", ";
    }
    std::cout << std::endl;
    for(int i = 0; i < numRanges; i++){
        srm_probMatrix[1][i] = (double) srm_countMatrix[1][i] / (double) totalCounts;
        std::cout << srm_probMatrix[1][i] << ", ";
    }
    std::cout << std::endl;

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


