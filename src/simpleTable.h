//
//  simpleNtable.h
//  StochHMM
//
//  Created by Paul Lott on 11/20/12.
//  Copyright (c) 2012 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#ifndef StochHMM_simpleNtable_h
#define StochHMM_simpleNtable_h
#include <vector>
#include <algorithm>
#include <iostream>
#include "sequences.h"

namespace StochHMM{
	
	
	
	
	template <typename T>
	class simpleTable
	{
	private:
		std::vector<T>* elements; // Points to all the actual elements
		size_t  num_of_elements; // Total number of array elements
		size_t  num_of_alphabets;
		size_t  num_of_dimensions;
		
		std::vector<size_t> order;
		std::vector<size_t> alpha;
		
		std::vector<std::pair<size_t,size_t>> upstream;
		
		std::vector<size_t>  dimensions; // Sizes of the N dimensions
		std::vector<size_t>  subarray_length; // Dimensions of subarrays
	public:
		
		simpleTable<T>(std::vector<size_t>& order, std::vector<size_t>& alpha, std::vector<size_t> tracks, T init): elements(NULL), num_of_elements(1){
			num_of_alphabets = alpha.size();
			
			if (num_of_alphabets != order.size()){
				std::cerr << "Differing number of orders and alphabets handed to SimpleArray" << std::endl;
				exit(2);
			}
			
			for(size_t i=0;i<num_of_alphabets;++i){
				num_of_elements*= pow(alpha[i],order[i]+1);
				for(size_t j=0;j<=order[i];++j){
					dimensions.push_back(alpha[i]);
					
					std::pair<size_t,size_t> temp(tracks[i],order[i]-j);
					upstream.push_back(temp);
				}
			}
			
			num_of_dimensions = dimensions.size();
			
			for (size_t i=0;i<num_of_dimensions;++i){
				size_t temp = 1;
				for (size_t j=i+1;j<num_of_dimensions;++j){
					temp*=dimensions[j];
				}
				subarray_length.push_back(temp);
			}
			
			elements = new (std::nothrow) std::vector<T>(num_of_elements, init);
			
			if (elements == NULL){
				std::cerr << "Unable to allocate SimpleArray" <<std::endl;
			}
			
			
			return;
		};
		
		
		~simpleTable<T>() { delete elements; }
		
		
		T get(sequences& seqs, size_t position){
			size_t index(0);
			
			for(size_t i=0; i < num_of_dimensions; i++){
				index += seqs[upstream[i].first][position-upstream[i].second] * subarray_length[i];
			}
			return (*elements)[index];
		}
		
		const T get(sequences& seqs, size_t position) const{
			size_t index(0);
			for(size_t i=0; i < num_of_dimensions; i++){
				index += seqs[upstream[i].first][position-upstream[i].second] * dimensions[i];
			}
			
			return elements[index];
		}
		
		void assign(T& value, size_t index){
			if (index > num_of_elements){
				std::cerr << "Adding element beyond range of array" << std::endl;
			}
			(*elements)[index] = value;
		};
		
		void assign(T& value, std::vector<uint8_t>& seq){
			size_t index(0);
			
			for(size_t i=0; i < num_of_dimensions; i++){
				index += seq[i] * dimensions[i];
			}
			
			(*elements)[index]=value;
		}
		
		size_t size()       const { return num_of_elements; }
		bool empty()           const { return num_of_elements==0; }
		
		// Return the size of each dimension
		size_t size(unsigned int Dim) const{ return dimensions[Dim-1];}
		
		// Delete all array elements
		void clear();
	};
	
	
	
	

}



#endif
