//emm.cpp
 //Copyright (c) 2007-2012 Paul C Lott 
 //University of California, Davis
 //Genome and Biomedical Sciences Facility
 //UC Davis Genome Center
 //Ian Korf Lab
 //Website: www.korflab.ucdavis.edu
 //Email: lottpaul@gmail.com
 //
 //Permission is hereby granted, free of charge, to any person obtaining a copy of
 //this software and associated documentation files (the "Software"), to deal in
 //the Software without restriction, including without limitation the rights to
 //use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 //the Software, and to permit persons to whom the Software is furnished to do so,
 //subject to the following conditions:
 //
 //The above copyright notice and this permission notice shall be included in all
 //copies or substantial portions of the Software.
 //
 //THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 //IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 //FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 //COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 //IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 //CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include "emm.h"
namespace StochHMM{

    //! Create an emission
    emm::emm(){
        real_number			= false;
        complement			= false;
        continuous			= false;
		multi_continuous	= false;
		realTrack			= NULL;
        
        function			= false;
        lexFunc				= NULL;
		
		pdf					= NULL;
		dist_parameters		= NULL;
		
		multiPdf			= NULL;
		number_of_tracks	= 0;
		trcks				= NULL;
		pass_values			= NULL;
		track_indices		= NULL;
		
        tagFunc				= NULL;
    }

        
    //! Destroy an emission    
    emm::~emm(){
        delete lexFunc;
        delete tagFunc;
        function =	false;

        lexFunc	=	NULL;
        tagFunc	=	NULL;
    }
    
    //!Parse an emission from text model file
    //!\param txt  String representation of emission
    //!\param trks Tracks used by the model
    //!\param wts Weights used by the model
    //!\param funcs State functions used by the model
    bool emm::parse(std::string& txt,tracks& trks, weights* wts, StateFuncs* funcs){
        if (!_processTags(txt,trks, wts, funcs)){
            return false;
        }
        
        stringList ln;
        ln.splitString(txt,"\n");
        size_t idx;
        if (ln.contains("EMISSION")){
            idx = ln.indexOf("EMISSION");
        }
        else{
            std::cerr << "Missing EMISSION tag from emission. Please check the formatting.   This is what was handed to the emission class:\n " <<  txt << std::endl;
            return false;
        }
        
        stringList line;
        line.splitString(ln[idx], "\t,: ");
        
        size_t typeBegin(0);
        
		//Determine Emission Type and set appropriate flags
        valueType  valtyp(PROBABILITY);
        if (line.contains("P(X)")){
            typeBegin = line.indexOf("P(X)");
            valtyp=PROBABILITY;
        }
        else if (line.contains("LOG")){
            typeBegin = line.indexOf("LOG");
            valtyp=LOG_PROB;
        }
        else if (line.contains("COUNTS")){
            typeBegin = line.indexOf("COUNTS");
            valtyp=COUNTS ;
        }
        else if (line.contains("REAL_NUMBER")){
            typeBegin = line.indexOf("REAL_NUMBER");
            real_number = true;
            if (line.contains("COMPLEMENT") || line.contains("1-P(X)")) {
                complement=true;
            }
        }
		else if (line.contains("MULTI_CONTINUOUS")){
			typeBegin = line.indexOf("MULTI_CONTINUOUS");
			multi_continuous=true;
			if (line.contains("COMPLEMENT") || line.contains("1-P(X)")) {
                complement=true;
            }
		}
		else if (line.contains("CONTINUOUS")){
			typeBegin = line.indexOf("CONTINUOUS");
			continuous=true;
			if (line.contains("COMPLEMENT") || line.contains("1-P(X)")) {
                complement=true;
            }
		}
		
        else if (line.contains("FUNCTION")){
            typeBegin = line.indexOf("FUNCTION");
            function=true;
        }
        else {
            std::string info = "Couldn't parse Value type in the Emission: " + txt  + " Please check the formatting.   The allowed types are: P(X), LOG, COUNTS, or REAL_NUMBER. \n";
            std::cerr << info << std::endl;
            
            //errorInfo(sCantParseLine, info.c_str());
        }
        
        
        //remaining tracks and Orders then set Track
        std::vector<track*> tempTracks;
        for(size_t i=1;i<typeBegin;i++){
            track* tk = trks.getTrack(line[i]);
            if (tk==NULL){
                std::cerr << "Emissions tried to add a track named: " << line[i] << " . However, there isn't a matching track in the model.  Please check to model formatting.\n";
                return false;
            }
            else{
                tempTracks.push_back(tk);
            }
        }
        
        //Real Number Emissions
        if (real_number){
            
            if (tempTracks.size()>1){
                std::cerr << "Multiple tracks listed under Real Track Emission Definition\n";
                return false;
            }
            
            realTrack = tempTracks[0];
            return true;
        }
		//Multivariate Continuous PDF emission
		else if (multi_continuous){
			if (tempTracks.size()==1){
				std::cerr << "Only a single track listed under MULTI_CONTINUOUS\n\
				Use CONTINUOUS instead of MULTI-CONTINUOUS\n";
				return false;
			}
			
			//Assign track information
			track_indices = new std::vector<size_t>;
			trcks = new std::vector<track*> (tempTracks);
			number_of_tracks = trcks->size();
			pass_values = new std::vector<double> (number_of_tracks,-INFINITY);
			for(size_t i = 0; i < number_of_tracks ; ++i){
				track_indices->push_back((*trcks)[i]->getIndex());
			}
			
			idx = ln.indexOf("PDF");
			line.splitString(ln[idx],"\t:, ");
			
			size_t function_idx = line.indexOf("PDF") + 1;
			multiPdfName = line[function_idx];
            multiPdf = funcs->getMultivariatePdfFunction(multiPdfName);
			
			
			size_t parameter_idx = line.indexOf("PARAMETERS");
			
			dist_parameters = new(std::nothrow) std::vector<double>;
			
			if(parameter_idx != SIZE_MAX){
				for(size_t i = parameter_idx+1 ; i< line.size() ; i++){
					double value;
					stringToDouble(line[i], value);
					dist_parameters->push_back(value);
				}
			}
			
			return true;
			
		}
		//U
		else if (continuous){
			
			if (tempTracks.size()>1){
                std::cerr << "Multiple tracks listed under CONTINUOUS Track Emission Definition\n\
				Must use MULTI_CONTINUOUS for multivariate emissions\n";
                return false;
            }
            
            realTrack = tempTracks[0];
			
			idx = ln.indexOf("PDF");
			line.splitString(ln[idx],"\t:, ");
			
			size_t function_idx = line.indexOf("PDF") + 1;
			pdfName = line[function_idx];
			
			size_t parameter_idx = line.indexOf("PARAMETERS");
			dist_parameters = new(std::nothrow) std::vector<double>;
			
			
			if(parameter_idx != SIZE_MAX){
				for(size_t i = parameter_idx+1 ; i< line.size() ; i++){
					double value;
					stringToDouble(line[i], value);
					dist_parameters->push_back(value);
				}
			}
			
            pdf = funcs->getPDFFunction(pdfName);
            
            return true;
		}
        else if (function){
			//Get function name
            std::string& functionName = line[typeBegin+1];
			
			//Set parameters for function
            lexFunc = new(std::nothrow) emissionFuncParam(functionName,funcs,tempTracks[0]);
            
            if (lexFunc==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
            
            return true;
        }
        else{  //Traditional Lexcical Emission
            
            
            if (ln.contains("ORDER")){
                idx=ln.indexOf("ORDER");
            }
            else{
                std::cerr << "Couldn't find ORDER in non-Real_Number emission.  Please check the formatting" << std::endl;
                return false;
                //errorInfo(sCantParseLine, "Couldn't find ORDER in non-Real_Number emission.  Please check the formatting\n");
            }
            
            std::vector<int> tempOrder;
            line.splitString(ln[idx],"\t:, ");
            
            size_t orderIdx = line.indexOf("ORDER");
            orderIdx++;
            
            size_t ambIdx; 
            bool containsAmbig = line.contains("AMBIGUOUS");
            
            if (containsAmbig){ambIdx=line.indexOf("AMBIGUOUS");}
            else{ ambIdx= line.size();}
            
            for(size_t i=orderIdx;i<ambIdx;i++){
                
                int tempValue;
                if (!stringToInt(line[i], tempValue)){
                    std::cerr << "Emission Order not numeric" << std::endl;
                    return false;
                }
                
                if (tempValue>32){
                    std::cerr << "Emission order is greater than 32.  Must be 32 or less" << std::endl;
                    return false;
                }
                
                tempOrder.push_back(tempValue);
            }
            
            if (tempOrder.size() == tempTracks.size()){
                for(size_t i=0;i<tempOrder.size();i++){
                    scores.addTrack(tempTracks[i], tempOrder[i]);
                }
            }
            else{
                std::cerr << "Different number of tracks and orders parsed in Emission: " << txt << " Check the formatting of the Emission" << std::endl;
                return false;
            }
            
            
            //Parse Ambiguous Tag Info
            if (containsAmbig){
                ambIdx++;
                if (line.size()<=ambIdx){
                    std::cerr << "No scoring type after AMBIGUOUS label\nAssuming AVG\n";
                    scores.setUnkScoreType(AVERAGE_SCORE);
                }
                else if (line[ambIdx].compare("AVG")==0){scores.setUnkScoreType(AVERAGE_SCORE);}
                else if (line[ambIdx].compare("MAX")==0){scores.setUnkScoreType(HIGHEST_SCORE);}
                else if (line[ambIdx].compare("MIN")==0){scores.setUnkScoreType(LOWEST_SCORE);}
				
				//Constant values assigned for Ambiguous characters
				//Can either be passed as P(X) or LOG
                else if (line[ambIdx].compare("P(X)")==0)
                {
                    scores.setUnkScoreType(DEFINED_SCORE);
                    
                    ambIdx++;
                    if (ambIdx>=line.size()){
                        std::cerr << "Missing Ambiguous Value" << std::endl;
                        
                        return false;
                    }
                                    
                    double tempValue;
                    if (!stringToDouble(line[ambIdx], tempValue)){
                        std::cerr << "Ambiguous Value couldn't be parsed: "<< line[ambIdx] << std::endl;
                        return false;
                    }
                    
                    scores.setUnkScore(log(tempValue));
                }
                else if (line[ambIdx].compare("LOG")==0){
                    scores.setUnkScoreType(DEFINED_SCORE);
                    ambIdx++;
                    
                    if (ambIdx>=line.size()){
                        std::cerr << "Missing Ambiguous Value" << std::endl;
                        
                        return false;
                    }
                    
                    double tempValue;
                    if (!stringToDouble(line[ambIdx], tempValue)){
                        std::cerr << "Ambiguous Value couldn't be parsed: "<< line[ambIdx] << std::endl;
                        return false;
                    }
                    
                    scores.setUnkScore(tempValue);
                }
            }
        
            //Get Emission Tables
            size_t expectedColumns(1);
            size_t expectedRows(1);
            for(size_t i = 0; i<scores.getNumberOfAlphabets(); i++){
                expectedColumns*=scores.getAlphaSize(i);
                expectedRows*=POWER[scores.getOrder(i)][scores.getAlphaSize(i)-1];
            }
            
            std::vector<std::vector<double> >* log_prob = scores.getLogProbabilityTable();
            std::vector<std::vector<double> >* prob = scores.getProbabilityTable();
            std::vector<std::vector<double> >* counts = scores.getCountsTable();
            
            
            for (size_t iter = 2; iter< ln.size();iter++){
                
                
                //If it's the first line check for a '#' indicating that the column header is present
                if (iter==2 && ln[iter][0]=='@'){
                    continue;
                }
                
                line.splitString(ln[iter],"\t ");
                
                //Check for Row header
                if (line[0][0]=='@'){
                    line.pop_ith(0);
                }
                
                
                std::vector<double> temp = line.toVecDouble();
                if (temp.size() != expectedColumns){
                    std::string info = "The following line couldn't be parsed into the required number of columns.   Expected Columns: " + int_to_string(expectedColumns) + "\n The line appears as: "  + ln[iter] ;
                    
                    std::cerr << info << std::endl;
                    return false;
                    //errorInfo(sCantParseLine, info.c_str());
                }
                else{
                    if (valtyp == PROBABILITY){
                        prob->push_back(temp);
                        logVector(temp);
                        log_prob->push_back(temp);
                    }
                    else if (valtyp == LOG_PROB){
                        log_prob->push_back(temp);
                        expVector(temp);
                        prob->push_back(temp);
                    }
                    else if (valtyp == COUNTS){
                        counts->push_back(temp);
                        probVector(temp);
                        prob->push_back(temp);
                        logVector(temp);
                        log_prob->push_back(temp);
                    }
                }
            }
            
            if (log_prob->size() != expectedRows){
                std::cerr << " The Emission table doesn't contain enough rows.  Expected Rows: " << expectedRows << " \n Please check the Emission Table and formatting for " <<  txt << std::endl;
                return false;
            }
			
			scores.initialize_emission_table();
			
        }
        
        return true;
    }
	
	//!Parse an emission from text
    //!\param txt  String representation of emission
    //!\param trks Tracks used by the model
    //!\param wts Weights used by the model
    //!\param funcs State functions used by the model
    bool emm::parse(std::string& txt,track* trk){
        
        stringList ln;
        ln.splitString(txt,"\n");
        size_t idx;
        if (ln.contains("EMISSION")){
            idx = ln.indexOf("EMISSION");
        }
        else{
            std::cerr << "Missing EMISSION tag from emission. Please check the formatting.   This is what was handed to the emission class:\n " <<  txt << std::endl;
            return false;
        }
        
        stringList line;
        line.splitString(ln[idx], "\t,: ");
		//size_t typeBegin(0);
                
        valueType  valtyp(PROBABILITY);
        if (line.contains("P(X)")){
            //typeBegin = line.indexOf("P(X)");
            valtyp=PROBABILITY;
        }
        else if (line.contains("LOG")){
            //typeBegin = line.indexOf("LOG");
            valtyp=LOG_PROB;
        }
        else if (line.contains("COUNTS")){
            //typeBegin = line.indexOf("COUNTS");
            valtyp=COUNTS ;
        }
		else {
            std::string info = "Couldn't parse Value type in the Emission: " + txt  + " Please check the formatting.   The allowed types are: P(X), LOG, COUNTS, or REAL_NUMBER. \n";
            std::cerr << info << std::endl;
            
            //errorInfo(sCantParseLine, info.c_str());
        }
        
        
        //remaining tracks and Orders then set Track
        std::vector<track*> temp_tracks;
		temp_tracks.push_back(trk);
        
		
		if (ln.contains("ORDER")){
			idx=ln.indexOf("ORDER");
		}
		else{
			std::cerr << "Couldn't find ORDER in non-Real_Number emission.  Please check the formatting" << std::endl;
			return false;
			//errorInfo(sCantParseLine, "Couldn't find ORDER in non-Real_Number emission.  Please check the formatting\n");
		}
		
		std::vector<int> temp_order;
		line.splitString(ln[idx],"\t:,");
		
		size_t orderIdx = line.indexOf("ORDER");
		orderIdx++;
		
		size_t ambIdx;
		bool containsAmbig = line.contains("AMBIGUOUS");
		
		if (containsAmbig){ambIdx=line.indexOf("AMBIGUOUS");}
		else{ ambIdx= line.size();}
		
		for(size_t i=orderIdx;i<ambIdx;i++){
			
			int temp_value;
			if (!stringToInt(line[i], temp_value)){
				std::cerr << "Emission Order not numeric" << std::endl;
				return false;
			}
			
			if (temp_value>32){
				std::cerr << "Emission order is greater than 32.  Must be 32 or less" << std::endl;
				return false;
			}
			
			temp_order.push_back(temp_value);
		}
		
		if (temp_order.size() == temp_tracks.size()){
			for(size_t i=0;i<temp_order.size();i++){
				scores.addTrack(temp_tracks[i], temp_order[i]);
			}
		}
		else{
			std::cerr << "Different number of tracks and orders parsed in Emission: " << txt << " Check the formatting of the Emission" << std::endl;
			return false;
		}
		
		
		//Parse Ambiguous Tag Info
		if (containsAmbig){
			ambIdx++;
			if (line.size()<=ambIdx){
				std::cerr << "No scoring type after AMBIGUOUS label\nAssuming AVG\n";
				scores.setUnkScoreType(AVERAGE_SCORE);
			}
			else if (line[ambIdx].compare("AVG")==0){scores.setUnkScoreType(AVERAGE_SCORE);}
			else if (line[ambIdx].compare("MAX")==0){scores.setUnkScoreType(HIGHEST_SCORE);}
			else if (line[ambIdx].compare("MIN")==0){scores.setUnkScoreType(LOWEST_SCORE);}
			else if (line[ambIdx].compare("P(X)")==0)
			{
				scores.setUnkScoreType(DEFINED_SCORE);
				
				ambIdx++;
				if (ambIdx>=line.size()){
					std::cerr << "Missing Ambiguous Value" << std::endl;
					
					return false;
				}
				
				double tempValue;
				if (!stringToDouble(line[ambIdx], tempValue)){
					std::cerr << "Ambiguous Value couldn't be parsed: "<< line[ambIdx] << std::endl;
					return false;
				}
				
				scores.setUnkScore(log(tempValue));
			}
			else if (line[ambIdx].compare("LOG")==0){
				scores.setUnkScoreType(DEFINED_SCORE);
				ambIdx++;
				
				if (ambIdx>=line.size()){
					std::cerr << "Missing Ambiguous Value" << std::endl;
					
					return false;
				}
				
				double tempValue;
				if (!stringToDouble(line[ambIdx], tempValue)){
					std::cerr << "Ambiguous Value couldn't be parsed: "<< line[ambIdx] << std::endl;
					return false;
				}
				
				scores.setUnkScore(tempValue);
			}
		}
		
		//Get Tables
		size_t expectedColumns(1);
		size_t expectedRows(1);
		for(size_t i = 0; i<scores.getNumberOfAlphabets(); i++){
			expectedColumns*=scores.getAlphaSize(i);
			expectedRows*=POWER[scores.getOrder(i)][scores.getAlphaSize(i)-1];
		}
		
		std::vector<std::vector<double> >* log_prob = scores.getLogProbabilityTable();
		std::vector<std::vector<double> >* prob = scores.getProbabilityTable();
		std::vector<std::vector<double> >* counts = scores.getCountsTable();
		
		
		for (size_t iter = 2; iter< ln.size();iter++){
			
			
			//If it's the first line check for a '#' indicating that the column header is present
			if (iter==2 && ln[iter][0]=='@'){
				continue;
			}
			
			line.splitString(ln[iter],"\t ");
			
			//Check for Row header
			if (line[0][0]=='@'){
				line.pop_ith(0);
			}
			
			
			std::vector<double> temp = line.toVecDouble();
			if (temp.size() != expectedColumns){
				std::string info = "The following line couldn't be parsed into the required number of columns.   Expected Columns: " + int_to_string(expectedColumns) + "\n The line appears as: "  + ln[iter] ;
				
				std::cerr << info << std::endl;
				return false;
				//errorInfo(sCantParseLine, info.c_str());
			}
			else{
				if (valtyp == PROBABILITY){
					prob->push_back(temp);
					logVector(temp);
					log_prob->push_back(temp);
				}
				else if (valtyp == LOG_PROB){
					log_prob->push_back(temp);
					expVector(temp);
					prob->push_back(temp);
				}
				else if (valtyp == COUNTS){
					counts->push_back(temp);
					probVector(temp);
					prob->push_back(temp);
					logVector(temp);
					log_prob->push_back(temp);
				}
			}
		}
		
		if (log_prob->size() != expectedRows){
			std::cerr << " The Emission table doesn't contain enough rows.  Expected Rows: " << expectedRows << " \n Please check the Emission Table and formatting for " <<  txt << std::endl;
			return false;
		}
		
		scores.initialize_emission_table();
		
		if (tagFunc != NULL){
			std::cerr << "Not NULL" << std::endl;
		}
        return true;
    }

	
	//!Parses the Emission Function Tag information from the text model definition
    bool emm::_processTags(std::string& txt, tracks& trks,weights* wts, StateFuncs* funcs){
        stringList lst = extractTag(txt);
                
        if (lst.size() == 0){
            return true;
        }
        
        tagFunc=new(std::nothrow) emissionFuncParam();
        
        
        if (tagFunc==NULL){
            std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
            exit(1);
        }
        
        
        if (!tagFunc->parse(lst, trks, wts,funcs)){
            std::cerr << "Couldn't parse Emission Tag: " << lst.stringify() << std::endl;
            return false;
        }
        std::string trackName = tagFunc->getTrackName();
        
        
        //Don't need if we check in the parsing of emissionFuncParam....
        if (!trks.isTrackDefined(trackName)){
            std::cerr << "No Track defined with name:\t" << trackName << "\nEmission Tag:\n" << txt << std::endl;
            return false;
        }
        
        
        return true;
    }


    
    //!Calculate the emission value given a position in the sequence
    //!If emission is a real number it will return the value from the real number track
    //!If emission is a sequence then it will get the value and return it
    //!\param seqs  Sequences to use
    //!\param iter Position within the sequences
    //!\return double log(prob) value of emission
    double emm::get_emission(sequences& seqs,size_t pos){
        double final_emission(-INFINITY);
        
        if (real_number){
            
            final_emission=seqs.realValue(realTrack->getIndex(),pos);
            if (complement){
                final_emission=log(1-exp(final_emission));
            }
        }
        else if (function){
            final_emission=lexFunc->evaluate(seqs, pos);
        }
		else if (multi_continuous){
			//Get all values from the tracks
			for (size_t i = 0; i < number_of_tracks ; ++i){
				(*pass_values)[i] = seqs.realValue((*track_indices)[i], pos);
			}
			
			final_emission = (*multiPdf)(pass_values, dist_parameters);
			
			if (complement){
				final_emission=log(1-exp(final_emission));
			}
		}
		else if (continuous){
			
			final_emission = (*pdf)(seqs.realValue(realTrack->getIndex(),pos),dist_parameters);
			
			if (complement){
				final_emission=log(1-exp(final_emission));
			}
		}
        else{
            final_emission=scores.getValue(seqs, pos);
        }
        
        
        if (tagFunc!=NULL){
            final_emission+=tagFunc->evaluate(seqs, pos);
        }
        
        return final_emission;
    }
	
	
	//!Calculate the emission value given a position in the sequence
    //!If emission is a real number it will return the value from the real number track
    //!If emission is a sequence then it will get the value and return it
    //!\param seqs  Sequences to use
    //!\param iter Position within the sequences
    //!\return double log(prob) value of emission
    double emm::get_emission(sequence& seq,size_t pos){
        double final_emission;
        
        if (real_number){
            
            final_emission=seq.realValue(pos);
            if (complement){
                final_emission=log(1-exp(final_emission));
            }
        }
        else if (function){
            final_emission=lexFunc->evaluate(seq, pos);
        }
		else if (multi_continuous){
			//Get all values from the tracks
			for (size_t i = 0; i < number_of_tracks ; ++i){
				(*pass_values)[i] = seq.realValue(pos);
			}
			
			final_emission = (*multiPdf)(pass_values, dist_parameters);
			
			if (complement){
				final_emission=log(1-exp(final_emission));
			}
		}
		else if (continuous){
			
			final_emission = (*pdf)(seq.realValue(pos),dist_parameters);
			
			if (complement){
				final_emission=log(1-exp(final_emission));
			}
		}
        else{
            final_emission=scores.getValue(seq, pos);
        }
        
        
        if (tagFunc!=NULL){
            final_emission+=tagFunc->evaluate(seq, pos);
        }
        
        return final_emission;
    }

    
    
    //!Is the emission from a real Number track
    //!\return true if track is real number track
    //!\return false if track is alphanumerical track
    bool emm::isReal(){
        if (real_number && realTrack->getAlphaType()==REAL){
            return true;
        }
        else{
            return false;
        }
    }

        
    //Get the string representation of the emission
    //\return std::string 
    std::string emm::stringify(){
        std::string emissionString("EMISSION:\t");
        
        if (real_number){
            
            emissionString+=realTrack->getName();
            emissionString+=":\t";
            
            emissionString+="REAL_NUMBER";
            if (complement){
                emissionString+=":\tCOMPLEMENT\t";
            }
            else{
                emissionString+="\t";
            }
            
            
            if (tagFunc){
                emissionString+=tagFunc->stringify();
            }
            emissionString+="\n";
        }
		//Univariate Continuous PDF emission 
		else if (continuous){
			emissionString+=realTrack->getName();
			emissionString+=":\tCONTINUOUS";
			
			if (tagFunc){
                emissionString+=tagFunc->stringify();
            }
            emissionString+="\n\t";
			
			emissionString+="PDF:\t";
			emissionString+=pdfName + "\tPARAMETERS:\t";
			emissionString+=join(*dist_parameters, ',');
			emissionString+= "\n";
			
		}
		else if (multi_continuous){
			for (size_t i=0; i < number_of_tracks; i++) {
				if (i>0){
					emissionString+=",";
				}
				emissionString+=(*trcks)[i]->getName();
			}
			emissionString+=":\tMULTI_CONTINUOUS";
			
			if (tagFunc){
                emissionString+=tagFunc->stringify();
            }
            emissionString+="\n\t";
			
			emissionString+="PDF:\t";
			emissionString+=pdfName + "\tPARAMETERS:\t";
			emissionString+=join(*dist_parameters, ',');
			emissionString+= "\n";
		}
        else if (function){
            emissionString+=lexFunc->getTrack()->getName();
            emissionString+=":\tFUNCTION:\t";
            emissionString+=lexFunc->getName();
            emissionString+="\t";
            
            if (tagFunc){
                emissionString += tagFunc->stringify();
            }
            
            emissionString+="\n";
        }
        else{
            
            for(size_t i=0;i<scores.trackSize();i++){
                if (i>0){
                    emissionString+=",";
                }
                emissionString+=scores.getTrack(i)->getName();
            }
            
            emissionString+=":\t";
            
            emissionString+="LOG";
            
            if (tagFunc){
                
                emissionString += "\t";
                
                emissionString += tagFunc->stringify();
            }
            
            emissionString+="\n\tORDER:\t";
            
            for(size_t i=0;i<scores.trackSize();i++){
                if (i>0){
                    emissionString+=",";
                }
                emissionString+=int_to_string(scores.getOrder(i));
            }
            
            unknownCharScoringType  ambTemp = scores.getAmbScoringType();
            
            if (ambTemp!=NO_SCORE){
                emissionString+="\tAMBIGUOUS:\t";
                emissionString+=(ambTemp==HIGHEST_SCORE)? "MAX":
                                (ambTemp==LOWEST_SCORE)? "MIN":
                                (ambTemp==AVERAGE_SCORE)? "AVG": "LOG:" + double_to_string(scores.getAmbDefinedScore());
            }
            emissionString+="\n";

            emissionString+=scores.stringify();
			//scores.stringifyAmbig();
        }
        emissionString+="\n";
        
        return emissionString;
    }

}
