//
//  Lexical.cpp
//  StochHMM
//
//  Created by Paul Lott on 4/2/12.
//  Copyright (c) 2012 University of California, Davis. All rights reserved.
//

#include "lexicalTable.h"



namespace StochHMM{
    
    lexicalTable::lexicalTable(){
        max_order=0;
        
        logProb=NULL;
        counts = NULL;
        prob = NULL;
		log_emission = NULL;
        
        
        unknownScoreType=NO_SCORE;
        unknownDefinedScore=-INFINITY;
        
        return;
    }
    
    lexicalTable::~lexicalTable(){
        delete logProb;
        delete prob;
        delete counts;
		delete log_emission;
        
        logProb=NULL;
        prob=NULL;
        counts=NULL;
		log_emission = NULL;
    }
    
    void lexicalTable::createTable(int rows, int columns, int pseudocount, valueType typ){
        if (typ==COUNTS){
            if (counts!=NULL){
                delete counts;
            }
            std::vector<double> temp_columns(columns,pseudocount);
            counts=new(std::nothrow) std::vector<std::vector<double> > (rows,temp_columns);
            
            if (counts==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
        }
        return;
    }
    
    void lexicalTable::print() {
        std::cout << stringify() << std::endl;
    }
    
    
    std::vector<std::vector<double> >* lexicalTable::getCountsTable(){
        if (counts==NULL){
            counts = new(std::nothrow) std::vector<std::vector<double> >;
            
            if (counts==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
			
        }
        
        return counts;
    }
    
    
    std::vector<std::vector<double> >* lexicalTable::getProbabilityTable(){
        if (prob==NULL){
            prob = new(std::nothrow) std::vector<std::vector<double> >;
            
            if (prob==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
			
        }
        
        return prob;
    }
    
    
    std::vector<std::vector<double> >* lexicalTable::getLogProbabilityTable(){
        if (logProb==NULL){
            logProb = new(std::nothrow) std::vector<std::vector<double> >;
            
            if (logProb==NULL){
                std::cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << std::endl;
                exit(1);
            }
        }
        
        return logProb;
    }
    
    //Return emission probability of sequences
    double lexicalTable::getValue(sequences& seqs, size_t pos){
        
		if (max_order>pos){
			return getReducedOrder(seqs, pos);
		}
		
		size_t index(seqs[subarray_sequence[0]][pos - subarray_position[0]] * subarray_value[0]);
		
		for(size_t i=1;i<dimensions;i++){
			index += seqs[subarray_sequence[i]][pos - subarray_position[i]] * subarray_value[i];
		}
		
		if (index > array_size){
			std::cerr << "Index is out of range of lookup table in lexicalTable" << std::endl;
			exit(2);
		}
		
		return (*log_emission)[index];
    }
	
	//Return emission probability of sequences
    double lexicalTable::getValue(sequence& seq, size_t pos){
        
		if (max_order>pos){
			return getReducedOrder(seq, pos);
		}
		
		size_t index(seq[pos - subarray_position[0]] * subarray_value[0]);
		
		for(size_t i=1;i<dimensions;i++){
			index += seq[pos - subarray_position[i]] * subarray_value[i];
		}
		
		if (index > array_size){
			std::cerr << "Index is out of range of lookup table in lexicalTable" << std::endl;
			exit(2);
		}
		
		return (*log_emission)[index];
    }
    
    
    //!Add a track to an emission
    //!\param trk Pointer to track
    //!\param orderValue order of emission from track
    void lexicalTable::addTrack(track* trk,int orderValue){
        trcks.push_back(trk);
        alphabets.push_back(trk->getAlphaSize());
        order.push_back(orderValue);
        if (orderValue>max_order){
            max_order=orderValue;
        }
        
    }
	
    //!Set the type of counts in the emission 2D table provided by the user
    //!\param temp vector of vectors of doubles
    //!\param emmType Type of value (COUNTS, PROBABILITY, LOG_PROB)
    void lexicalTable::assignTable(std::vector<std::vector<double > >* temp, valueType emmType){
        if (emmType==COUNTS){
            counts=temp;
        }
        else if (emmType == PROBABILITY){
            prob=temp;
        }
        else if (emmType == LOG_PROB){
            logProb=temp;
			initialize_emission_table();
        }
    }
	
    
    std::string lexicalTable::stringify(){
        std::string tbl("");
        size_t tracks_size = trcks.size();
        
        if(tracks_size==0){
            std::cerr << "Can't print out table without track and order being set for lexicalTable\n";
            exit(1);
        }
        
        //Output Column Headers
        size_t columns(1);
        std::vector<size_t> alphaSizes;
        //alphaSizes.push_back(0);
        
        for(size_t i = 0;i<trcks.size();i++){
            size_t alphaSz = trcks[i]->getAlphaSize();
            columns*=alphaSz;
            alphaSizes.push_back(alphaSz);
        }
        
        
        reverse(alphaSizes.begin(),alphaSizes.end());
        
        std::string colHeader("@");
        
        for(size_t i = 0;i<columns;i++){
            size_t indexValue = i;
            for(size_t tr=0;tr<trcks.size();tr++){
                
                if (tr>0){
                    colHeader+= "|";
                }
                
                size_t val(0);
                if (tr<trcks.size()-1){
                    val= floor(indexValue/alphaSizes[tr]);
                    indexValue-=val*alphaSizes[tr];
                }
                else{
                    val = indexValue;
                }
                
                colHeader+=trcks[tr]->convertIndexToWord(val, 1);
            }
            colHeader+="\t";
        }
        
        tbl+=colHeader + "\n";
        
        std::vector<std::vector<double> >* temp;
        
        if (logProb!=NULL){
            temp=logProb;
        }
        else if (prob!=NULL){
            temp=prob;
        }
        else if (counts!=NULL){
            temp=counts;
        }
        else{
            std::cerr << "No table is defined\n";
            return "";
        }
        
        //TODO: Fix row header for other orders
        bool rowHeader = (temp->size()>1) ? true : false;
        for(size_t i=0;i<temp->size();i++){
            std::string header("");
            
            if (rowHeader){
                size_t indexValue = i;
                
                for(size_t tr=0;tr<trcks.size();tr++){
                    
                    if (tr>0 && order[tr]>0){
                        header+= "|";
                    }
                    
                    
                    size_t val(0);
                    if (tr<trcks.size()-1){
                        double pwr = POWER[order[tr+1]][trcks[tr+1]->getAlphaSize()-1];
                        val= floor(indexValue/pwr);
                        indexValue-=val*pwr;
                    }
                    else{
                        val = indexValue;
                    }
                    
                    header+=trcks[tr]->convertIndexToWord(val, order[tr]);
                }
                tbl+="@" + header + "\t";
                
            }
            
            for(size_t j=0;j<(*temp)[i].size();j++){
                if (j>0){
                    tbl+="\t";
                }
                tbl+=double_to_string((*temp)[i][j]);
            }
            tbl+="\n";
        }
        return tbl;
    }
	
	
	void lexicalTable::init_table_dimension_values(){
		number_of_tracks = trcks.size();
		y_dim = sumVector(order);
		
		
		//Calculate subarray dimensions for logProb and new log_emission table (includes ambiguous character)
		x_subarray = new size_t[number_of_tracks];
		y_subarray = new size_t[y_dim];
		
		//Calculate Old subarray x_dimension values
		for(size_t i=0;i<number_of_tracks;++i){
			size_t value(1);
			for(size_t j=i+1;j<number_of_tracks;++j){
				value*=alphabets[j];
			}
			x_subarray[i]=value;
		}
		
		//Calcuate Old subarray y_dimension values
		std::vector<size_t> index(order[0],alphabets[0]);
		for(size_t i=1;i<order.size();i++){
			for(size_t j=0;j<order[i];j++){
				index.push_back(alphabets[i]);
			}
		}
		
		for (size_t i=0;i<y_dim;i++){
			size_t val = 1;
			for(size_t j=i+1;j<y_dim;j++){
				val*=index[j];
			}
			y_subarray[i]=val;
		}
		
		return;
	}
	
	//TODO:  NEED TO COMPLETE THIS FUNCTION
//	size_t lexicalTable::convertIndex(size_t old_row, size_t old_column){
//		//Decompose previous value from indices to digital letter value
//		
//	}
	
	
	void lexicalTable::init_array_dimension_values(){
		dimensions = y_dim + number_of_tracks;
		subarray_sequence.resize(dimensions);
		subarray_value.resize(dimensions);
		subarray_value.resize(dimensions);
		
		decompose_values.resize(dimensions);
		decompose_sequence.resize(dimensions);
		
		//Calculate total size of emission table needed with ambiguous characters
		array_size = 1;
		std::vector<size_t> complete_alphabet_size;
		size_t current_dim(0);
		for(size_t i=0;i<number_of_tracks;i++){
			
			//Get alphabet size and store it
			size_t alpha_size = trcks[i]->getTotalAlphabetSize();
			complete_alphabet_size.push_back(alpha_size);
			array_size*=integerPower(alpha_size, (size_t) order[i]+1);
			
			for(size_t j=0;j<=order[i];++j){
				subarray_sequence[current_dim]=i;
				current_dim++;
			}
		}
		
		//Calculate the Sequence positions
		for (size_t i=0;i<number_of_tracks;i++){
			for (size_t j=0;j<=order[i];++j){
				subarray_position.push_back(order[i]-j);
			}
		}
		

		//Compute the decomposing values
		//Used to convert from index to word
		std::vector<size_t> decompose_index;
		for(size_t i=0;i<number_of_tracks;++i){
			for(size_t j=0;j<order[i];++j){
				decompose_index.push_back(complete_alphabet_size[i]);
			}
		}
		for(size_t i=0;i<number_of_tracks;++i){
			decompose_index.push_back(complete_alphabet_size[i]);
		}
		
		//Calculate subarray values
		for(size_t i=0;i<dimensions;++i){
			decompose_values[i]=1;
			for(size_t j=i+1;j<dimensions;++j){
				decompose_values[i]*=decompose_index[j];
			}
		}
		
		
		//Rearrange decompose values for subarray values;
		//Final values  in subarray_value will correspond to sequence AAA(A)B(B).
		//Where A is 3rd order and B is 1st order;
		size_t array_index(0);
		size_t index(0);
		for(size_t i=0;i<number_of_tracks;++i){
			for(size_t j=0;j<order[i];++j){
				subarray_value[array_index] = decompose_values[index];
				array_index++;
				index++;
			}
			subarray_value[array_index] = decompose_values[dimensions-number_of_tracks-i];
			array_index++;
		}
		
		
		//Finalize decompose_sequence
		std::vector<size_t> temp;
		for(size_t i=0;i<number_of_tracks;++i){
			for(size_t j=0;j<order[i];++j){
				temp.push_back(i);
			}
		}
		for(size_t i=0;i<number_of_tracks;++i){
			temp.push_back(i);
		}
		
		for(size_t i=0;i<dimensions;++i){
			decompose_sequence[i]=temp[i];
		}
		
		return;
	}
	
	
	//Transfer values from 2d table to array and also compute the ambiguous score
	void lexicalTable::transferValues(std::vector<bool>& transferred){
		
		//Transfer unambiguous scores
		for(size_t row=0;row<logProb->size();++row){
			for(size_t column=0;column<(*logProb)[row].size();++column){
				std::vector<uint8_t> alphabet;
				decompose(row, column, alphabet);
				size_t index = calculateArrayIndex(alphabet);
				(*log_emission)[index] = (*logProb)[row][column];
				transferred[index] = true;
			}
		}
		
		//Compute all ambiguous scores
		for(size_t i=0;i<transferred.size();i++){
			if (transferred[i]){
				continue;
			}
			
			if (unknownScoreType == DEFINED_SCORE){
				(*log_emission)[i] = unknownDefinedScore;
				continue;
			}
			else if (unknownScoreType == NO_SCORE){
				continue;
			}
			
			//Get the letters for index
			std::vector<uint8_t> letters;
			decompose(i,letters);
			
			//Expand the unambiguous words and get all the possible values
			//std::vector<double> expanded;
			//expand_ambiguous(letters,expanded);
			
			(*log_emission)[i] = getAmbiguousScore(letters);
			
//			//Assign the values accordint the Score Type
//			if (unknownScoreType == AVERAGE_SCORE){
//				(*log_emission)[i] = avgLogVector(expanded);
//			}
//			else if (unknownScoreType == LOWEST_SCORE){
//				(*log_emission)[i] = minVector(expanded);
//			}
//			else if (unknownScoreType == HIGHEST_SCORE){
//				(*log_emission)[i] = maxVector(expanded);
//			}
		}
		
//		for (size_t i=0;i<log_emission->size();++i){
//			std::vector<uint8_t> letters;
//			decompose(i,letters);
//			for (size_t j = 0; j< letters.size();j++){
//				std::cout << (int) letters[j];
//			}
//			std::cout << "\t" ;
//			std::cout << (*log_emission)[i] << std::endl;
//		}
		return;
	}
	
	//Calculate lower order emission from the current table values
	//Given order and position/sequence
	//Calculate the values using Index and [all alphabets] for higher orders
	double lexicalTable::getReducedOrder(sequences& seqs, size_t position){
		Index indices;
		for(size_t i=0;i<dimensions;i++){
			Index subtotal;
			size_t seq = subarray_sequence[i];
			size_t pos = subarray_position[i];
			
			if (subarray_position[i] > position){
				subtotal.setAmbiguous(trcks[seq]->getUnambiguousSet());
			}
			else if (seqs[seq][position - pos] > max_unambiguous[seq]){
				subtotal.setAmbiguous(trcks[seq]->getAmbiguousSet(seqs[seq][position-pos]));
			}
			else{
				subtotal+=seqs[seq][position-pos];
			}
			
			subtotal *= subarray_value[i];
			indices  += subtotal;
		}
		
		if (unknownScoreType == AVERAGE_SCORE || unknownScoreType == NO_SCORE){
			double temp(0);
			for(size_t i=0;i<indices.size();i++){
				temp+=exp((*log_emission)[indices[i]]);
			}
			temp /= indices.size();
			temp = log(temp);
			return temp;
		}
		else if (unknownScoreType == LOWEST_SCORE){
			double temp(INFINITY);
			for(size_t i=0;i<indices.size();i++){
				double val = (*log_emission)[indices[i]];
				if (val < temp){
					temp = val;
				}
			}
			return temp;
		}
		else if (unknownScoreType == HIGHEST_SCORE){
			double temp(-INFINITY);
			for(size_t i=0;i<indices.size();i++){
				double val = (*log_emission)[indices[i]];
				if (val > temp){
					temp = val;
				}
			}
			return temp;
		}
		return -INFINITY;
	}
	
	
	//Calculate lower order emission from the current table values
	//Given order and position/sequence
	//Calculate the values using Index and [all alphabets] for higher orders
	double lexicalTable::getReducedOrder(sequence& seq, size_t position){
		Index indices;
		for(size_t i=0;i<dimensions;i++){
			Index subtotal;
			size_t sq = subarray_sequence[i];
			size_t pos = subarray_position[i];
			
			if (subarray_position[i] > position){
				subtotal.setAmbiguous(trcks[sq]->getUnambiguousSet());
			}
			else if (seq[position - pos] > max_unambiguous[sq]){
				subtotal.setAmbiguous(trcks[sq]->getAmbiguousSet(seq[position-pos]));
			}
			else{
				subtotal+=seq[position-pos];
			}
			
			subtotal *= subarray_value[i];
			indices  += subtotal;
		}
		
		if (unknownScoreType == AVERAGE_SCORE || unknownScoreType == NO_SCORE){
			double temp(0);
			for(size_t i=0;i<indices.size();i++){
				temp+=exp((*log_emission)[indices[i]]);
			}
			temp /= indices.size();
			temp = log(temp);
			return temp;
		}
		else if (unknownScoreType == LOWEST_SCORE){
			double temp(INFINITY);
			for(size_t i=0;i<indices.size();i++){
				double val = (*log_emission)[indices[i]];
				if (val < temp){
					temp = val;
				}
			}
			return temp;
		}
		else if (unknownScoreType == HIGHEST_SCORE){
			double temp(-INFINITY);
			for(size_t i=0;i<indices.size();i++){
				double val = (*log_emission)[indices[i]];
				if (val > temp){
					temp = val;
				}
			}
			return temp;
		}
		return -INFINITY;
	}

	
	double lexicalTable::getAmbiguousScore(std::vector<uint8_t>& letters){
		Index indices;
		for(size_t i=0;i<dimensions;++i){
			Index subtotal;
			if (letters[i]>max_unambiguous[decompose_sequence[i]]){
				subtotal.setAmbiguous(trcks[decompose_sequence[i]]->getAmbiguousSet(letters[i]));
			}
			else{
				 subtotal+= letters[i];
			}
			
			subtotal *= decompose_values[i];
			indices += subtotal;
		}
		
		if (unknownScoreType == AVERAGE_SCORE){
			double temp(0);
			for(size_t i=0;i<indices.size();i++){
				temp+=exp((*log_emission)[indices[i]]);
			}
			temp /= indices.size();
			temp = log(temp);
						
			return temp;
		}
		else if (unknownScoreType == LOWEST_SCORE){
			double temp(INFINITY);
			for(size_t i=0;i<indices.size();i++){
				double val = (*log_emission)[indices[i]];
				if (val < temp){
					temp = val;
				}
			}
			return temp;
		}
		else if (unknownScoreType == HIGHEST_SCORE){
			double temp(-INFINITY);
			for(size_t i=0;i<indices.size();i++){
				double val = (*log_emission)[indices[i]];
				if (val > temp){
					temp = val;
				}
			}
			return temp;
		}
		
		return -INFINITY;
	}
	
	
	//Use index instead much faster than expanding the letters and then reverse computing 
//	void lexicalTable::expand_ambiguous(std::vector<uint8_t>& letters, std::vector<double>& expanded){
//		
//		std::vector<std::vector<uint8_t> >* temp_words = new std::vector<std::vector<uint8_t> >;
//		temp_words->push_back(letters);
//		for(size_t i=0;i<dimensions;i++){
//			temp_words = expand_ambiguous(temp_words, i);
//		}
//		
//		for(size_t i=0;i<temp_words->size();i++){
//			size_t index = calculateIndexFromDecomposed((*temp_words)[i]);
//			expanded.push_back((*log_emission)[index]);
//		}
//		
//		return;
//	}
//	
//	std::vector<std::vector<uint8_t> >* lexicalTable::expand_ambiguous(std::vector<std::vector<uint8_t> >* words, size_t letter){
//		
//		std::vector<std::vector<uint8_t> >* temp_words= new std::vector<std::vector<uint8_t> >;
//		for(size_t i=0; i < words->size(); ++i){
//			if ((*words)[i][letter] > max_unambiguous[decompose_sequence[letter]]){
//				std::vector<uint8_t>& set = trcks[decompose_sequence[letter]]->getAmbiguousSet((*words)[i][letter]);
//				for(size_t j=0;j<set.size();++j){
//					(*words)[i][letter] = set[j];
//					temp_words->push_back((*words)[i]);
//				}
//			}
//			else{
//				temp_words->push_back((*words)[i]);
//			}
//		}
//		
//		delete words;
//		words = NULL;
//		return temp_words;
//	}
	
	size_t lexicalTable::calculateIndexFromDecomposed(std::vector<uint8_t>& word){
		size_t index(0);
		for(size_t i=0;i<dimensions;++i){
			index += decompose_values[i] * word[i];
		}
		return index;
	}
	
	size_t lexicalTable::calculateArrayIndex(std::vector<uint8_t>& kmer){
		size_t index(0);
		for(size_t i=0;i<dimensions;++i){
			index += subarray_value[i] * kmer[i];
		}
		return index;
	}
	
	void lexicalTable::decompose(size_t row, size_t column, std::vector<uint8_t>& letters){
		//Decompose the row into the preceding letters
		for(size_t i=0;i<y_dim;++i){
			size_t val = floor(row/y_subarray[i]);
			row -= val*y_subarray[i];
			letters.push_back(val);
		}
		
		//Decompose the column into the emitted letters
		for(size_t i=0;i<number_of_tracks;++i){
			size_t val = floor(column/x_subarray[i]);
			column-=val*x_subarray[i];
			letters.push_back(val);
		}
		
		return;
	}
	
	
	void lexicalTable::decompose(size_t index, std::vector<uint8_t>& letters){
		
		//Decompose the index into the emitted letters
		for(size_t i=0;i<dimensions;++i){
			size_t val = floor(index/subarray_value[i]);
			index-=val*subarray_value[i];
			letters.push_back(val);
		}
		
		return;
	}
	
	//Todo
	//Convert table to simpleNtable compatible format
	//with ambiguous characters
	void lexicalTable::initialize_emission_table(){
		if (logProb == NULL){
			std::cerr << "Cannot initialize emission table until after the tables have been assigned";
			exit(2);
		}
		init_table_dimension_values();
		init_array_dimension_values();
		
		for(size_t i = 0; i < number_of_tracks ; ++i){
			max_unambiguous.push_back(trcks[i]->getMaxUnambiguous());
		}
		
		if(unknownDefinedScore == DEFINED_SCORE){
			log_emission = new std::vector<double> (array_size,unknownDefinedScore);
		}
		else{
			log_emission = new std::vector<double> (array_size,-INFINITY);
		}
		
		//Transfer values to emission_table
		std::vector<bool> transferred (array_size,false);
		transferValues(transferred);
	}
	
    
    
}