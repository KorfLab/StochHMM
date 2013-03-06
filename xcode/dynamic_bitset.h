//
//  dynamic_bitset.h
//  dynamic_bitset
//
//  Created by Paul Lott on 3/1/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#ifndef __dynamic_bitset__dynamic_bitset__
#define __dynamic_bitset__dynamic_bitset__

#include <iostream>
#include <vector>
#include <math.h>
#include <stdint.h>
#include "bitwise_ops.h"

//!\class dynamic_bitset
//! A dynamic bitset class
//! Allows the bitset to be set dynamically and manipulated in size
//! dynamic_bitset uses 32 bit unsigned integers as the underlying
//! datatype.
//! Need to implement shifting to right and left
//! Note: This dynamic_bitset is read from left to right
//! This class can be greatly improved by using SSE optimization
//! For portability and time issues: SSE optimization haven't been implemented yet.
class dynamic_bitset{
public:
	dynamic_bitset():current_size(0),buffer(0),num_ints(0){}
	dynamic_bitset(size_t);
	dynamic_bitset(const dynamic_bitset&);
	dynamic_bitset(const std::string&);
	//	dynamic_bitset(std::vector<unsigned char>&);
	//	dynamic_bitset(std::vector<uint16_t>&);
	//	dynamic_bitset(std::vector<uint32_t>&);
	//	dynamic_bitset(std::vector<uint64_t>&);
	//	dynamic_bitset(size_t val, uint16_t*);
	//	dynamic_bitset(size_t val, uint32_t*);
	//	dynamic_bitset(size_t val, uint64_t*);
	//	dynamic_bitset(size_t val, unsigned char*);
	
	void resize(size_t);
	void reserve(size_t);
	inline bool empty(){return (current_size ==0) ? true : false;}
	void clear();
	
	bool operator[](size_t) const;
	dynamic_bitset  operator&	(const dynamic_bitset& rhs);
	dynamic_bitset  operator|	(const dynamic_bitset& rhs);
	dynamic_bitset  operator^	(const dynamic_bitset& rhs);
	dynamic_bitset  operator~	() const;
	
	dynamic_bitset& operator=	(const dynamic_bitset& rhs);
	dynamic_bitset& operator&=	(const dynamic_bitset& rhs);
	dynamic_bitset& operator|=	(const dynamic_bitset& rhs);
	dynamic_bitset& operator^=	(const dynamic_bitset& rhs);
	dynamic_bitset& operator-=	(const dynamic_bitset& rhs);
	dynamic_bitset& operator<<=	(size_t n);
	dynamic_bitset& operator>>=	(size_t n);
	dynamic_bitset  operator<<	(size_t n);
	dynamic_bitset  operator>>	(size_t n);
	
	void insert(size_t pos, size_t n);
	void erase(size_t pos);
	
	dynamic_bitset get_intersection(const dynamic_bitset& rhs);
	dynamic_bitset get_union(const dynamic_bitset& rhs);
	dynamic_bitset get_unique(const dynamic_bitset& rhs);
	
	bool operator== (const dynamic_bitset& rhs) const;
	bool operator!= (const dynamic_bitset& rhs) const;
	
	friend std::ostream& operator<< (std::ostream& , const dynamic_bitset&);
	
	inline size_t size(){return current_size;}
	
	bool at(size_t pos)const;
	bool test(size_t pos);
	
	void set(size_t pos);
	void set(size_t pos, bool value);
	
	void unset(size_t pos);
	
	void flip();
	void flip(size_t pos);
	
	void reset();
	void reset(size_t pos);
	
	bool any();
	bool none();
	
	void push_back(bool);
	size_t count();
	size_t count_before(size_t);
	
	size_t find_first() const;
	size_t find_first(size_t pos) const;
	
	size_t find_last() const;
	size_t find_last(size_t pos) const;
	
	bool parity();
	
	std::string stringify() const;
	std::string stringify_all() const;
	
	bool intersects(const dynamic_bitset& )const;
	
private:
	size_t buffer;
	size_t current_size;
	size_t num_ints;
	//uint32_t* array;
	std::vector<uint32_t> array;
	
};

#endif /* defined(__dynamic_bitset__dynamic_bitset__) */
