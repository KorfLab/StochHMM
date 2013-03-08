//
//  sparseArray.h
//  DynamicBitset
//
//  Created by Paul Lott on 2/28/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#ifndef __DynamicBitset__sparseArray__
#define __DynamicBitset__sparseArray__

#include <iostream>
#include <vector>
#include "dynamic_bitset.h"

namespace StochHMM {

	/*!SparseArray is a template class that stores the values in a sparse array
	 * The values in the array are indexed by a dynamic bitset.  For example, to get
	 the value at iterator 5, we need to convert that iterator to the iterator of the
	 std::vector array.  This is essentially the number of values that are stored 
	 before the value at the 5th position.
	 
	 SparseArray uses std::vector.   This means that it will incure an increased cost
	 of inserting and deleting values.
	 
	 */
	template <typename T>
	class sparseArray{
	public:
		sparseArray();
		sparseArray(size_t sz):flag(sz){};
		sparseArray(const sparseArray& val):flag(val.flag), array(val.array){};
		
		void clear();
		void reserve(size_t);
		void resize(size_t);
		
		//!Returns the number of elements stored in the array
		//TODO: which is faster (Population count on bitset or get size from vector)
		inline size_t elements(){return flag.count();};
		
		//!Returns the size of the sparseArray
		inline size_t size(){return flag.size();}
		
		
		//Reference to Elements in Array
		T& operator[] (size_t);
		const T& operator[] (size_t) const;
		T& at(size_t pos){};
		T& back();
		const T& back()const;
		
		//Comparison Operators
		bool operator==	(const sparseArray& rhs);
		bool operator!=	(const sparseArray& rhs);
		
		//Assignment(Copy) Operator
		sparseArray& operator= (const sparseArray& rhs);
		
		//!Return whether an element is defined at position
		//!\returns true if element is defined at position
		//!\returns false if element is not defined at position
		inline bool defined(size_t pos){return flag[pos];}
		
		//Adding Elements to sparseArray
		void insert(size_t pos, T& val);
		void push_back(T&);
		
		//Removing Elements from sparseArray
		void erase(size_t pos);
		void pop_back(size_t pos);

	private:
		dynamic_bitset flag;
		std::vector<T> array;
	};

	//!Removes all the elements from the array
	template <typename T>
	void sparseArray<T>::clear() {
		flag.clear();
		array.clear();
		return;
	}

	//!Reserved memory for at least n elements in the array;
	template <typename T>
	void sparseArray<T>::reserve(size_t n) {
		flag.reserve(n);
		array.reserve(n);
		return;
	}

	//!Resize the sparseArray to n;
	template <typename T>
	void sparseArray<T>::resize(size_t n) {
		//Increasing size only affect increase in flag size
		if (n >= flag.size()){
			flag.resize(n);
		}
		// Decrease size- need to also delete items from the array
		else{
			flag.resize(n);
			//# of values to keep is equal to bits set in flag
			size_t sz = flag.count();
			
			//Resize array to values needed
			array.resize(sz);
		}
		return;
	}

	//!Pushes the element on to the array
	template <typename T>
	void sparseArray<T>::push_back(T& val){
		flag.push_back(true);
		array.push_back(val);
		return;
	}

	//!Insert an element at the position in the array
	//!\param pos Position within the array to insert value
	//!\param reference to value to insert
	template <typename T>
	void sparseArray<T>::insert(size_t pos,T& val){
		flag.insert(pos, 1);
		flag.set(pos);
		size_t array_iter = flag.count_before(pos);
		array.insert(array.begin()+array_iter,val);
		return;
	}

	//!Returns reference to object at position
	//!If the position is beyond the length, the sparse array will get extended.
	//!If it is defined it will return a reference to the position
	//!If it is not defined and within the current length of the sparseArray and
	//! after all the currently set values , it will get added to the end of the array
	//!If it is not defined and within the current lenght of the sparseArray and
	//! before other currently set values, then it will get inserted into the array
	template <typename T>
	T& sparseArray<T>::operator[](size_t pos){
		//Increase size of flags if necessary
		if (flag.size() <= pos){
			flag.resize(pos);
		}
		
		//if called on defined data then return reference
		if (defined(pos)){
			return array[flag.count_before(pos)];
		}
		
		//if called on undefined data then that value will be added as set and returned
		//with default value of type
		
		//If position is beyond the last set position then we can simply add it to
		//the array
		if (flag.find_last() < pos){
			flag.set(pos);
			array.resize(array.size()+1);
			return array.back();
		}
		
		//Otherwise, we need to insert the position into the array
		//Uses default constructor of typename T
		flag.set(pos);
		size_t array_iter = flag.count_before(pos);
		array.insert(array.begin()+array_iter,T());
		
		return array[array_iter];
	}


	//!Returns reference to object at position
	//!If the position is beyond the length, the sparse array will get extended.
	//!If it is defined it will return a reference to the position
	//!If it is not defined and within the current length of the sparseArray and
	//! after all the currently set values , it will get added to the end of the array
	//!If it is not defined and within the current lenght of the sparseArray and
	//! before other currently set values, then it will get inserted into the array
	template <typename T>
	const T& sparseArray<T>::operator[](size_t pos)const {
		
		//Increase size of flags if necessary
		if (flag.size() <= pos){
			flag.resize(pos);
		}
		
		//if called on defined data then return reference
		if (defined(pos)){
			return array[flag.count_before(pos)];
		}
		
		//if called on undefined data then that value will be added as set and returned
		//with default value of type
		
		//If position is beyond the last set position then we can simply add it to
		//the array
		if (flag.find_last() < pos){
			flag.set(pos);
			array.resize(array.size()+1);
			return array.back();
		}
		
		//Otherwise, we need to insert the position into the array
		//Uses default constructor of typename T
		flag.set(pos);
		size_t array_iter = flag.count_before(pos);
		array.insert(array.begin()+array_iter,T());
		
		return array[array_iter];
	}


}


#endif /* defined(__DynamicBitset__sparseArray__) */
