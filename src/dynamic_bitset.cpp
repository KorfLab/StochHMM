//
//  dynamic_bitset.cpp
//  dynamic_bitset
//
//  Created by Paul Lott on 3/1/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "dynamic_bitset.h"

namespace StochHMM {
		
	dynamic_bitset::dynamic_bitset(size_t sz): current_size(sz){
		num_ints = ((sz-1)/32)+1;
		buffer = (num_ints*32) - sz;
		array.resize(num_ints);
		//array = new (std::nothrow) uint32_t[num_ints]{0};
		return;
	}

	dynamic_bitset::dynamic_bitset(const dynamic_bitset& rhs){
		num_ints = rhs.num_ints;
		buffer	= rhs.buffer;
		current_size = rhs.current_size;
		array = rhs.array;
		return;
	}

	dynamic_bitset::dynamic_bitset(const std::string& str){
		size_t str_size = str.size();
		num_ints = ((str_size-1)/32)+1;
		buffer = (num_ints*32) - str_size;
		current_size=str_size;
		array.resize(num_ints);
		
		for(size_t i = 0; i < str_size; ++i){
			if (str[(str_size-1) - i ] == '1'){
				set(i);
				continue;
			}
			else if (str[(str_size-1) - i ] == '0'){
				continue;
			}
			else{
				std::cerr << "Invalid argument: " << str[(str_size-1) - i ] << " Only 0's and 1's \
				allowed in string to initialize dynamic_bitset" << std::endl;
				exit(2);
			}
		}
		return;
	}

	dynamic_bitset::dynamic_bitset(const std::vector<bool>& vec){
		size_t vec_size = vec.size();
		num_ints = ((vec_size-1)/32)+1;
		buffer	= (num_ints*32) - vec_size;
		current_size = vec_size;
		array.resize(num_ints);
		for (size_t i = 0; i < vec_size; ++i){
			if (vec[i]){
				set(i);
			}
		}
		return;
	}

	//dynamic_bitset::dynamic_bitset(unsigned char ch){
	//	num_ints = 1;
	//	buffer = 24;
	//	current_size = 8;
	//	array.resize(num_ints);
	//	array[0] = (ch << 24);
	//	return;
	//}

	//dynamic_bitset::dynamic_bitset(uint16_t val){
	//	num_ints = 1;
	//	buffer = 16;
	//	current_size = 16;
	//	array.resize(num_ints);
	//	array[0] = (val << 16);
	//	return;
	//}
	//
	//dynamic_bitset::dynamic_bitset(uint32_t val){
	//	num_ints = 1;
	//	buffer = 0;
	//	current_size = 32;
	//	array.resize(num_ints);
	//	array[0] = val;
	//	return;
	//}
	//
	//dynamic_bitset::dynamic_bitset(uint64_t val){
	//	num_ints = 2;
	//	buffer = 0;
	//	current_size = 64;
	//	array.resize(num_ints);
	//	array[0] = (uint32_t) ((val & 0xFFFFFFFF00000000) >> 32);
	//	array[1] = (uint32_t) (val & 0xFFFFFFFF);
	//	return;
	//}


	//Need to implement getting smaller
	//! Resize the bitset.
	//! If size is smaller than bitset, those values will get deleted
	//! If size is larger than bitset, it will get extended.
	void dynamic_bitset::resize(size_t sz){
		size_t new_num_ints = ((sz-1)/32)+1;
		
		if (num_ints == new_num_ints){
			current_size = sz;
		}
		else{
			array.resize(new_num_ints);
			current_size = sz;
			num_ints = new_num_ints;
		}
		
		//Set new buffer size
		buffer = (num_ints*32)-sz;
		
		//Clear bit in buffer
		//Clear those bits that are outside of the current scope
		uint32_t mask = (1 << 32 - buffer) - 1;
		array[num_ints-1] &= mask;
		return;
	}


	//!Reserves the set amount of memory
	//!Doesn't change any of the bitset size. Only reserve a set amount of memory.
	void dynamic_bitset::reserve(size_t sz){
		if (sz <= current_size + buffer){
			return;
		}
		else{
			size_t new_num_ints = ((sz-1)/32)+1;
			array.reserve(new_num_ints);
		}
		
		return;
	}

	//!Clears the complete bitset
	//!Size is reduced to zero
	void dynamic_bitset::clear(){
		buffer=0;
		current_size=0;
		num_ints=0;
		array.clear();
		return;
	}


	//!Bitwise - AND operator
	//!\returns a dynamic_bitset that is *this AND rhs
	dynamic_bitset dynamic_bitset::operator&	(const dynamic_bitset& rhs){
		dynamic_bitset output(*this);
		output&=rhs;
		return output;
	}

	//!Bitwise - OR operator
	//!\returns a dynamic_bitset that is *this OR rhs
	dynamic_bitset dynamic_bitset::operator|	(const dynamic_bitset& rhs){
		dynamic_bitset output(*this);
		output|=rhs;
		return output;
	}

	//!Bitwise - XOR operator
	//!\returns a dynamic_bitset that is *this XOR rhs
	dynamic_bitset dynamic_bitset::operator^	(const dynamic_bitset& rhs){
		dynamic_bitset output(*this);
		output^=rhs;
		return output;
	}

	//!Bitwise ~ Complement
	//!\returns dynamic_bitset that is complement of *this
	dynamic_bitset dynamic_bitset::operator~() const{
		dynamic_bitset output(*this);
		for(size_t i = 0; i < num_ints; ++i){
			output.array[i] = ~output.array[i];
		}
		
		//Clear those bits that are outside of the current scope
		uint32_t mask = (1 << 32-output.buffer) - 1;
		output.array[num_ints-1] &= mask;
		
		return output;
	}

	//!Overloaded assignment operator makes *this a copy of rhs
	//!\param rhs dynamic_bitset to copy
	//!\return *this;
	dynamic_bitset& dynamic_bitset::operator=	(const dynamic_bitset& rhs){
		num_ints = rhs.num_ints;
		buffer	= rhs.buffer;
		current_size = rhs.current_size;
		array = rhs.array;
		return *this;
	}

	//! Bitwise AND all the bits with all the bits in the rhs
	dynamic_bitset& dynamic_bitset::operator&=	(const dynamic_bitset& rhs){
		
		if (rhs.current_size != this->current_size){
			std::cerr << __FUNCTION__ << " called on dynamic_bitsets of different sizes." <<std::endl;
			exit(2);
		}
		
		for(size_t i = 0; i < num_ints; ++i){
			array[i] &= rhs.array[i];
		}
		
		//Clear those bits that are outside of the current scope
		uint32_t mask = (1 << 32-buffer) - 1;
		array[num_ints-1] &= mask;
		
		return *this;
	}

	//! Bitwise OR all the bits with all the bits in the rhs
	dynamic_bitset& dynamic_bitset::operator|=	(const dynamic_bitset& rhs){
		
		if (rhs.current_size != this->current_size){
			std::cerr << __FUNCTION__ << " called on dynamic_bitsets of different sizes." <<std::endl;
			exit(2);
		}
		
		
		for(size_t i = 0; i < num_ints; ++i){
			array[i] |= rhs.array[i];
		}
		
		
		//Clear those bits that are outside of the current scope
		uint32_t mask = (1 << 32-buffer) - 1;
		array[num_ints-1] &= mask;
		
		return *this;
	}

	//! Bitwise XOR all the bits with all the bits in the rhs
	dynamic_bitset& dynamic_bitset::operator^=	(const dynamic_bitset& rhs){
		
		if (rhs.current_size != this->current_size){
			std::cerr << __FUNCTION__ << " called on dynamic_bitsets of different sizes." <<std::endl;
			exit(2);
		}
		
		for(size_t i = 0; i < num_ints; ++i){
			array[i] ^= rhs.array[i];
		}
		
		
		//Clear those bits that are outside of the current scope
		uint32_t mask = (1 << 32-buffer) - 1;
		array[num_ints-1] &= mask;
		
		return *this;
	}


	//!Overloaded ostream for using dynamic_bitsets in ostream
	//!Usage: std::cout << bs << ....;
	std::ostream& operator<< (std::ostream& output , const dynamic_bitset& bs){
		output << bs.stringify();
		
		return output;
	}

	//!Checks to see if bitset is equal to *this
	//! \param Bitset to compare with *this
	//! \returns false if the bitsets are not equal to each other
	//! \returns true if the bitsets are equal
	bool dynamic_bitset::operator==(const dynamic_bitset& rhs) const{
		for(size_t i = 0; i < num_ints ;++i){
			if (array[i]!=rhs.array[i]){
				return false;
			}
		}
		return true;
	}

	//!Checks to see if bitset is not equal to *this
	//! \param Bitset to compare with *this
	//! \returns true if the bitsets are not equal to each other
	//! \returns false if the bitsets are equal
	bool dynamic_bitset::operator!=(const dynamic_bitset& rhs) const{
		for(size_t i = 0; i < num_ints ;++i){
			if (array[i]!=rhs.array[i]){
				return true;
			}
		}
		return false;
	}

	//!Shift the bitset to the left n positions
	//!Size of the bitset will not change.  If values are shifted off bitset then
	//!they will be lost.
	//!\param n Number of positions to shift the bits
	dynamic_bitset& dynamic_bitset::operator<<= (size_t n){
		//If shifting more than 32 then need to insert integer before
		size_t start(0);
		if (n>=32){
			size_t ints_to_add = n/32;
			n %=32;
			array.insert(array.begin(), ints_to_add, 0);
			array.resize(num_ints);
			if (n==0){
				return *this;
			}
		}
		
		if (n!=0){
			uint32_t shifted(lso(array[start],n));
			for(size_t i = start+1; i < num_ints ; ++i){
				shifted = lsoso(array[i], shifted, n);
			}
			
			//Clear those bits that are outside of the current scope (bits in buffer region)
			uint32_t mask = (1 << 32-buffer) - 1;
			array[num_ints-1] &= mask;
		}
		
		return *this;
	}

	dynamic_bitset dynamic_bitset::operator<< (size_t n){
		dynamic_bitset output(*this);
		output<<=n;
		return output;
	}


	//!Shift the bitset to the right n positions
	//!Size of the bitset will not change.  If values are shifted off bitset then
	//!they will be lost.
	//!\param n Number of positions to shift the bits
	dynamic_bitset& dynamic_bitset::operator>>= (size_t n){
		//If shifting more than 32 then need to insert integer before
		size_t start(num_ints-1);
		if (n>=32){
			size_t ints_to_delete = n/32;
			n %= 32;
			array.erase(array.begin(), array.begin()+ints_to_delete);
			start -= ints_to_delete;
			if (n == 0){
				array.resize(num_ints);
				return *this;
			}
		}
		
		if (n != 0){
			//Shift first integers and get value to shift to next
			uint32_t shifted(rso(array[start],n));
			
			//For all the integers before shift the value
			for(size_t i = start-1; i != SIZE_MAX ; --i){
				shifted = rsoso(array[i], shifted, n);
			}
			
			//Resesize the array
			array.resize(num_ints);
		}
		
		return *this;
	}

	dynamic_bitset dynamic_bitset::operator>> (size_t n){
		dynamic_bitset output(*this);
		output>>=n;
		return output;
	}

	//!Insert n bits (set to zero) at the position.  Bitsets at this position and
	//!greater will be shifted to the left.  This will increase the size of the bitset
	// \param pos	Position to insert the bits
	// \param n		Number of bits to insert
	void dynamic_bitset::insert(size_t pos, size_t n){
		if (pos > current_size){
			std::cerr << "dynamic_bitset::"<<__FUNCTION__ << " Error: Position for insert is beyond the range current size\n";
			exit(2);
		}
		
		size_t new_buffer = 32-(current_size + n)%32;
		size_t remainder = 32 - (pos%32);
		size_t new_num_ints = ((current_size-1)+n)/32+1;
		
		size_t first_affected = pos/32;
		//size_t last_affected = num_ints;
		
		size_t next_affected = (pos+n)/32;
		
		int64_t next_shift_size = ((pos+n)%32) - (pos%32);
		
		int64_t intervening_ints = int64_t (n-remainder) / 32;
		size_t relative_pos = pos%32;
		
		
		if (new_num_ints!= num_ints){
			
			if (intervening_ints != 0 && intervening_ints >= 0){
				array.insert(array.begin()+first_affected+1,intervening_ints,0);
				num_ints+= intervening_ints;
			}
			
			//Add additional on end if necessary
			if (new_num_ints != num_ints){
				array.resize(new_num_ints);
				num_ints= new_num_ints;
			}
		}
		
		uint32_t temp_mask = maskLeft(array[first_affected], relative_pos);
		array[first_affected] = maskRight(array[first_affected], relative_pos-1);
		
		//If shift is by factor of 32 then we just need to reassign values
		if (next_shift_size == 0){
			size_t ints_shifted = n/32;
			for (size_t i = num_ints-1; i > first_affected + intervening_ints + 1; --i){
				array[i] = array[i-(ints_shifted-intervening_ints)];
			}
			array[first_affected+ints_shifted] = temp_mask;
		}
		else if (next_shift_size < 0){ //Right shift after intervening sequence
			
			//Shift temp sequence right by necessary amount
			uint8_t left_shift = 32 - abs( (uint8_t) next_shift_size);
			temp_mask >>= abs( (uint8_t) next_shift_size);
			
			for(size_t i = next_affected; i < num_ints ; ++i){
				uint32_t tmp = array[i];
				array[i] <<= left_shift;
				array[i] |= temp_mask;
				temp_mask = tmp >> (32-left_shift);
				
				//temp_mask = lsoso(array[i], temp_mask, left_shift);
			}
			
		}
		else{//left_shift after intervening sequence
			
			//Move the values to new cells before shifting old values
			if (n>32){
				size_t ints_shifted = n/32;
				for (size_t i = num_ints-1; i > first_affected + intervening_ints + 1; --i){
					array[i] = array[i-(ints_shifted-intervening_ints)];
					array[i-(ints_shifted-intervening_ints)]=0;
				}
			}
			
			//Shift temp value
			uint32_t shifted = lso(temp_mask,next_shift_size);
			array[next_affected] |= temp_mask;
			
			//Shift the rest of the values
			for(size_t i = next_affected+1; i < num_ints ; ++i){
				uint32_t tmp = array[i];
				array[i] <<= next_shift_size;
				array[i] |= shifted;
				shifted = tmp >> (32-next_shift_size);
			}
		}
		buffer = new_buffer;
		current_size +=n;
		return;
	}


	//!Erase the bit at the position from the set.
	//!Causes the bits to shift right by one

	void dynamic_bitset::erase(size_t pos){
		
		//if the remove position is beyond the end the don't do anything.
		if (pos>=current_size){
			return;
		}
		
		//If size is equal to one then just assign zero and return
		if (current_size == 1){
			current_size =0;
			buffer = 32;
			array[0] = 0;
			return;
		}
		
		size_t new_buffer = 32-(current_size-1)%32;
		size_t new_num_ints = ((current_size-2))/32+1;
		size_t first_affected = pos/32;
		size_t relative_pos = pos%32;
		
		//Get the value shifted off the last int
		uint32_t shifted(0);
		
		//If not last
		if (first_affected != num_ints-1){
			shifted = rso(array[num_ints-1], 1);
			//Shift the rest of the integers up to first affected
			for (size_t i = num_ints-2; i != first_affected && i != SIZE_MAX; --i){
				uint32_t temp(array[i]);
				array[i] = (array[i]>>1) | shifted;
				shifted = temp << 31;
			}
		}
		
		
		//Get the left most part of the first affected
		uint32_t temp_mask = maskLeft(array[first_affected], relative_pos+1);
		//std::cout << intToBinString(temp_mask) << std::endl;
		
		//Get the right most part of the first affected and return it
		if ((relative_pos - 1) == SIZE_MAX){
			array[first_affected] = 0;
		}
		else{
			array[first_affected] = maskRight(array[first_affected], relative_pos-1);
		}
		
		
		//Shift the left part
		temp_mask >>=1;
		//std::cout << intToBinString(temp_mask) << std::endl;
		//Add the shifted bit from the previous
		temp_mask |= shifted;
		//std::cout << intToBinString(temp_mask) << std::endl;
		//Add pack to the first affected
		array[first_affected] |= temp_mask;
		
		//Resize the int if necessary
		if (new_num_ints != num_ints){
			array.resize(new_num_ints);
		}
		
		//Set new values of buffer and size;
		buffer = new_buffer;
		current_size--;
		
		return;
	}



	//!Return the intersection bitset with rhs bitset
	//! Same as bitwise AND
	//! \param dynamic_bitset&
	//! \return dynamic_bitset that is intersection of both
	dynamic_bitset dynamic_bitset::get_intersection(const dynamic_bitset& rhs){
		return *this & rhs;
	}

	//!Returns the union bitset
	//! Same as bitwise OR
	dynamic_bitset dynamic_bitset::get_union(const dynamic_bitset& rhs){
		return *this | rhs;
	}

	//!Return only the bits that are uniqe to bitset
	dynamic_bitset dynamic_bitset::get_unique(const dynamic_bitset& rhs){
		return *this &  ~rhs ;
	}




	//!Get state of bit at position
	//!No initial range check
	//!\returns true if set 1
	//!\returns false if not set 0
	bool dynamic_bitset::operator[](size_t pos)const{
		size_t intIter = pos/32;
		size_t bitIter = pos - (intIter*32);
		return array[intIter] & (1 << bitIter);
	}

	dynamic_bitset::bit_ref dynamic_bitset::operator[](size_t pos){
		size_t intIter = pos/32;
		size_t bitIter = pos - (intIter*32);
		return bit_ref(array[intIter],bitIter);
	}

	//!Returns state of bit at position
	//!Performs range check first
	bool dynamic_bitset::at(size_t pos) const{
		if (pos>=current_size){
			std::cerr << "Error: Out of Range.\nFunction: " <<  __FUNCTION__ << std::endl;
			exit(2);
		}
		
		size_t intIter = pos/32;
		size_t bitIter = pos - (intIter*32);
		return array[intIter] & (1 << bitIter);
	}

	//!Returns whether bit is set at the position
	//!Performs range check first
	//! \param pos Position within bitset
	//! \return true if bit is set
	//! \return false if bit is unset
	bool dynamic_bitset::test(size_t pos){
		if (pos>=current_size){
			std::cerr << "Error: Out of Range.\nFunction: " <<  __FUNCTION__ << std::endl;
			exit(2);
		}
		
		return at(pos);
	}

	//!Sets the bit at position with value
	//!Checks size before assignment
	//!\param position within bitset
	void dynamic_bitset::set(size_t pos){
		
		if (pos>=current_size){
			std::cerr << "Error: Out of Range.\nFunction: " <<  __FUNCTION__ << std::endl;
			exit(2);
		}
		
		size_t intIter = pos/32;
		size_t bitIter = pos - (intIter*32);
		array[intIter] |= (1 << bitIter);
		return;
	}

	//!Sets the bit at position with value
	//!Checks size before assignment
	//!\param position within bitset
	//!\param val Value to set bit
	void dynamic_bitset::set(size_t pos, bool val){
		if (pos>=current_size){
			std::cerr << "Error: Out of Range.\nFunction: " <<  __FUNCTION__ << std::endl;
			exit(2);
		}
		
		if (val){
			set(pos);
		}
		else{
			unset(pos);
		}
	}

	//!Unsets the bit at position
	//!Checks size before assignment
	void dynamic_bitset::unset(size_t pos){
		
		if (pos>=current_size){
			std::cerr << "Error: Out of Range.\nFunction: " <<  __FUNCTION__ << std::endl;
			exit(2);
		}
		
		size_t intIter = pos/32;
		size_t bitIter = pos - (intIter*32);
		array[intIter] &= ~(1 << bitIter);
		return;
	}

	//!Flips the whole bitset
	//!Complete bitset gets assigned the complement bitset
	void dynamic_bitset::flip(){
		*this = ~*this;
		return;
	}

	//!Flip the bit at position
	//! If bit is set 1, then it gets unset 0
	void dynamic_bitset::flip(size_t pos){
		if (test(pos)){
			unset(pos);
		}
		else{
			set(pos);
		}
		
		return;
	}

	//!Resets all the bits to zero
	void dynamic_bitset::reset(){
		array.assign(num_ints, 0);
	}

	//!Clears the bit as position
	void dynamic_bitset::reset(size_t pos){
		unset(pos);
		return;
	}

	//!Checks to see if any bits are set
	//!\returns true if any bits are set
	//!\returns false no bits are set
	bool dynamic_bitset::any(){
		for(size_t i = 0; i < num_ints; i++){
			if (array[i]!=0){
				return true;
			}
		}
		return false;
	}

	//!Checks to see if no bits are set
	//!\returns false if any bits are set
	//!\returns true if no bits are set
	bool dynamic_bitset::none(){
		for(size_t i = 0; i < num_ints; i++){
			if (array[i]!=0){
				return false;
			}
		}
		return true;
	}

	//!Pushes a bit onto the bitset
	//!Increases the bitset size by 1
	//!\param bool Value to set in bitset
	void dynamic_bitset::push_back(bool val){
		if (buffer==0){
			//Reserve more memory. Push uint32(0) onto array
			if (current_size != 0){
				reserve(current_size*2);
				array.push_back(0);
				buffer = 32;
				num_ints++;
			}
			else{
				reserve(32);
				array.push_back(0);
				buffer = 32;
				num_ints++;
			}
		}
		
		//Increases the size and decrease buffer;
		current_size++;
		buffer--;
		
		//Set value of position
		if (val){
			set(current_size-1);
		}
		else{
			unset(current_size-1);
		}
		
		return;
	}


	//!Find first set bit within bitset and return the position
	//!Searches from beginning to end.  (Left to Right)
	//!\return size_t position of first bit set
	//!\ref http://bits.stephan-brumme.com/lowestBitSet.html
	size_t dynamic_bitset::find_first() const{
		for(size_t i = 0; i < num_ints ; ++i){
			if (array[i] == 0){
				continue;
			}
			else{
				uint32_t val = array[i];
				
				static const size_t MultiplyDeBruijnBitPosition[32] =
				{
					// precomputed lookup table
					0,  1, 28,  2, 29, 14, 24,  3, 30, 22, 20, 15, 25, 17,  4,  8,
					31, 27, 13, 23, 21, 19, 16,  7, 26, 12, 18,  6, 11,  5, 10,  9
				};
				
				// leave only lowest bit
				val  &= -int(val);
				// DeBruijn constant
				val  *= 0x077CB531;
				// get upper 5 bits
				val >>= 27;
				// convert to actual position
				
				return i*32 + MultiplyDeBruijnBitPosition[val];
			}
		}
		return SIZE_MAX;
	}

	//!Find the first bit within the bitset starting at position pos
	//!Searches from beginning to end.  (Left to Right)
	//!\param pos Position to start the search
	//!\return Returns the position within the bitset with first set bit
	//\ref http://bits.stephan-brumme.com/lowestBitSet.html
	size_t dynamic_bitset::find_first(size_t pos) const {
		
		static const size_t MultiplyDeBruijnBitPosition[32] =
		{
			// precomputed lookup table
			0,  1, 28,  2, 29, 14, 24,  3, 30, 22, 20, 15, 25, 17,  4,  8,
			31, 27, 13, 23, 21, 19, 16,  7, 26, 12, 18,  6, 11,  5, 10,  9
		};
		
		if (pos > current_size){
			return SIZE_MAX;
		}
		
		//Deal with partial integer upto next integer boundary
		//Find Current interger
		size_t iter = pos/32;
		
		//Position within the integer to start
		size_t int_pos = pos%32;
		
		int32_t val = array[iter];
		
		//Mask bits before the position
		val = val & ~((1 << int_pos) - 1);
		
		if (val !=0){
			// leave only lowest bit
			val  &= -int(val);
			// DeBruijn constant
			val  *= 0x077CB531;
			// get upper 5 bits
			val >>= 27;
			// convert to actual position
			
			return iter * 32 +  MultiplyDeBruijnBitPosition[val];
		}
		
		//Deal with everything from integer boundary onward
		//http://bits.stephan-brumme.com/lowestBitSet.html
		for(size_t i = iter+1; i < num_ints ; ++i){
			if (array[i] == 0){
				continue;
			}
			else{
				uint32_t val = array[i];
				
				// leave only lowest bit
				val  &= -int(val);
				// DeBruijn constant
				val  *= 0x077CB531;
				// get upper 5 bits
				val >>= 27;
				// convert to actual position
				
				return i * 32 + MultiplyDeBruijnBitPosition[val];
			}
		}
		return SIZE_MAX;
	}


	//!Searches in reverse to find first set bit within bitset and return the position
	//!Searches from end to beginning. (Right to Left)
	//!\return size_t position of first bit set
	size_t dynamic_bitset::find_last() const{
		for(size_t i = num_ints-1; i != SIZE_MAX ; --i){
			if (array[i] == 0){
				continue;
			}
			else{
				//Adapted from http://stackoverflow.com/questions/671815/what-is-the-fastest-most-efficient-way-to-find-the-highest-set-bit-msb-in-an-i
				//Using binary search algorithm to search
				uint32_t val = array[i];
				const uint32_t mask[] = {
					0x00000FFFF,
					0x0000000FF,
					0x00000000F,
					0x000000003,
					0x000000001
				};
				uint32_t hi = 32;
				uint32_t lo = 0;
				
				if (val == 0)
					return 0;
				
				for (uint32_t j = 0; j < 5; j++) {
					int mi = lo + (hi - lo) / 2;
					
					if ((val >> mi) != 0)
						lo = mi;
					else if ((val & (mask[j] << lo)) != 0)
						hi = mi;
				}
				return lo + 32 * i;
			}
		}
		return SIZE_MAX;
	}


	size_t dynamic_bitset::find_last(size_t pos)const{
		if (pos > current_size){
			return SIZE_MAX;
		}
		
		//Deal with partial integer upto next integer boundary
		//Find Current interger
		size_t iter = pos/32;
		
		//Position within the integer to start
		size_t int_pos = pos%32;
		
		int32_t val = array[iter];
		//std::cout << intToBinString(val) << std::endl;
		
		//Mask bits before the position
		val = val & ((1 << int_pos+1) - 1);
		
		//std::cout << intToBinString(val) << std::endl;
		
		//If value isn't zero then we'll need to find which position is set within it
		if (val !=0) {
			//Adapted from http://stackoverflow.com/questions/671815/what-is-the-fastest-most-efficient-way-to-find-the-highest-set-bit-msb-in-an-i
			//Using binary search algorithm to search
			const uint32_t mask[] = {
				0x00000FFFF,
				0x0000000FF,
				0x00000000F,
				0x000000003,
				0x000000001
			};
			uint32_t hi = 32;
			uint32_t lo = 0;
			
			if (val == 0)
				return 0;
			
			for (uint32_t j = 0; j < 5; j++) {
				int mi = lo + (hi - lo) / 2;
				
				if ((val >> mi) != 0)
					lo = mi;
				else if ((val & (mask[j] << lo)) != 0)
					hi = mi;
			}
			return lo + 32 * iter;
		}
		
		for(size_t i = iter-1; i != SIZE_MAX ; --i){
			if (array[i] == 0){
				continue;
			}
			else{
				//Adapted from http://stackoverflow.com/questions/671815/what-is-the-fastest-most-efficient-way-to-find-the-highest-set-bit-msb-in-an-i
				//Using binary search algorithm to search
				uint32_t val = array[i];
				const uint32_t mask[] = {
					0x00000FFFF,
					0x0000000FF,
					0x00000000F,
					0x000000003,
					0x000000001
				};
				uint32_t hi = 32;
				uint32_t lo = 0;
				
				if (val == 0)
					return 0;
				
				for (uint32_t j = 0; j < 5; j++) {
					int mi = lo + (hi - lo) / 2;
					
					if ((val >> mi) != 0)
						lo = mi;
					else if ((val & (mask[j] << lo)) != 0)
						hi = mi;
				}
				return lo + 32 * i;
			}
		}
		
		return SIZE_MAX;
	}


	//!Returns whether there are an odd number of bits set
	//!\return true if number of 1's is odd
	//!\return false if number of 1' is even
	bool dynamic_bitset::parity(){
		bool val(false);
		for(size_t i = 0; i < num_ints; ++i){
			uint32_t temp = array[i];
			temp ^= temp >> 1;
			temp ^= temp >> 2;
			temp &= 0x011111111;
			temp *= 0x011111111;
			val ^= ((temp >> 28) & 1) != 0;
		}
		return val;
	}

	//!Returns the number of 1's in the bitset
	//!\ref bits.sephan-brumme countBits population
	//!\return size_t number of bits set in bitset
	size_t dynamic_bitset::count() const{
		size_t population(0);
		for (size_t i=0; i < num_ints; i++){
			uint32_t x = array[i];
			// count bits of each 2-bit chunk
			x  = x - ((x >> 1) & 0x55555555);
			// count bits of each 4-bit chunk
			x  = (x & 0x33333333) + ((x >> 2) & 0x33333333);
			// count bits of each 8-bit chunk
			x  = x + (x >> 4);
			// mask out junk
			x &= 0xF0F0F0F;
			// add all four 8-bit chunks
			population+= (x * 0x01010101) >> 24;
		}
		return population;
	}


	//!Returns the number of 1's in the bitset to the right of the current position
	//!\ref bits.sephan-brumme countBits population
	//!\return size_t number of bits set in bitset
	size_t dynamic_bitset::count_before(size_t pos) const{
		size_t current_int  = pos/32;
		size_t current_pos  = pos%32;
		
		uint32_t x(array[current_int]);
		
		//Mask all the bits to right (including pos);
		uint32_t mask = (1 << current_pos) - 1;
		x &= mask;
		
		//Count bits
		x = x - ((x >> 1) & 0x55555555);                    // reuse input as temporary
		x = (x & 0x33333333) + ((x >> 2) & 0x33333333);     // temp
		size_t population(((x + (x >> 4) & 0xF0F0F0F) * 0x1010101) >> 24); // count
		
		//Count the bit of integers to left
		for (size_t i=0; i < current_int; i++){
			x = array[i];
			// count bits of each 2-bit chunk
			x  = x - ((x >> 1) & 0x55555555);
			// count bits of each 4-bit chunk
			x  = (x & 0x33333333) + ((x >> 2) & 0x33333333);
			// count bits of each 8-bit chunk
			x  = x + (x >> 4);
			// mask out junk
			x &= 0xF0F0F0F;
			// add all four 8-bit chunks
			population+= (x * 0x01010101) >> 24;
		}
		return population;
	}


	//!Returns string representation of bitset
	//!\return std::string representation of bitset
	std::string dynamic_bitset::stringify() const{
		std::string output("");
		for (size_t i=current_size-1;i != SIZE_MAX; --i){
			
			
			bool test = at(i);
			if (test){
				output+="1";
			}
			else{
				output+="0";
			}
			
			//		if (i!= 0 && i%32 == 0){
			//			output+=",";
			//		}
			
		}
		return output;
	}

	//!Returns string representation of bitset includes buffered bits
	//\return std::string representation of bits including buffered bits
	std::string dynamic_bitset::stringify_all() const{
		std::string output("");
		for (size_t i=num_ints-1 ; i != SIZE_MAX; --i){
			output+= intToBinString(array[i]);
			//		if (i!=0){
			//			output+=",";
			//		}
		}
		
		return output;
	}


	dynamic_bitset::bit_ref& dynamic_bitset::bit_ref::operator=(bool value){
		//if true then set the bit at position
		if (value){
			val |= (1<<pos);
		}
		//else clear the bit at the position
		else{
			val &= ~(1<< pos);
		}
		return *this;
	}

	dynamic_bitset::bit_ref& dynamic_bitset::bit_ref::operator=(const bit_ref& rhs){
		//if set in rhs then set in this position
		if ((rhs.val >> pos) & 1){
			val |= (1 << pos);
		}
		//else clear the bit at this position
		else{
			val &= ~(1 << pos);
		}
		return *this;
	}

	dynamic_bitset::bit_ref& dynamic_bitset::bit_ref::operator|=(bool value){
		//if true then set bit
		if (value){
			val |= (1 << pos);
		}
		//Otherwise leave it alone
		return *this;
	}

	dynamic_bitset::bit_ref& dynamic_bitset::bit_ref::operator&=(bool value){
		//if value is false and val is true then clear val bit
		if (!value  && ((this->val >> pos) & 1)){
			val &= ~(1 << pos);
		}
		//otherwise leave unchanged
		return *this;
	}

	dynamic_bitset::bit_ref& dynamic_bitset::bit_ref::operator^=(bool value){
		//Flip the bit
		val^= (1 << pos);
		return *this;
	}

	dynamic_bitset::bit_ref& dynamic_bitset::bit_ref::flip(){
		//Flip the bit
		val^= (1 << pos);
		return *this;
	}

	dynamic_bitset::bit_ref& dynamic_bitset::bit_ref::reset(){
		//Clear the bit
		val &= ~(1 << pos);
		return *this;
	}

	bool dynamic_bitset::bit_ref::operator~()const{
		return !((this->val >> pos) & 1);
	}



}

