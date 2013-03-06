//
//  bitwise_ops.h
//  DynamicBitset
//
//  Created by Paul Lott on 3/1/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#ifndef __DynamicBitset__bitwise_ops__
#define __DynamicBitset__bitwise_ops__

#include <iostream>
#include <math.h>
#include <stdint.h>

std::string intToBinString(uint32_t val);
size_t msb(uint32_t val);
size_t lowestBitSet(uint32_t x);
int8_t popCount(uint32_t val);
uint32_t maskLeft(uint32_t val,uint8_t pos);
uint32_t maskRight(uint32_t val, uint8_t pos);

uint32_t rol(uint32_t& val, uint8_t shift_register);
uint32_t ror(uint32_t& val, uint8_t shift_register);

//!Rotate Left with carry over
template <typename T>
T rol(T val, uint8_t shift_register) {
    return (val << shift_register) | (val >> ((sizeof(T)*CHAR_BIT) - shift_register));
};


//!Rotate Right with carry over
template <typename T>
T ror(T val, uint8_t shift_register) {
    return (val >> shift_register) | (val << ((sizeof(T)*CHAR_BIT) - shift_register));
}

//!Left Shift with shift off
//!Shifts the value to the left by n bits.
//! The bits that would get shifted off value are returned in lsb
template <typename T>
T lso(T& val, uint8_t n) {
	T temp(val);
	if (n>=(sizeof(T)*CHAR_BIT)){
		val = 0;
	}
	else{
		val <<=n;
	}
	n%=32;
    return temp >> ((sizeof(T)*CHAR_BIT) - n);
};

//!Right Shift with shift off
//!Shifts the value to the right by n bits.
//! The bits that would get shifted off value are returned in msb
template <typename T>
T rso(T& val, uint8_t n) {
	T temp(val);
	if (n>=(sizeof(T)*CHAR_BIT)){
		val = 0;
	}
	else{
		val >>=n;
	}
	n%=32;
    return temp << ((sizeof(T)*CHAR_BIT) - n);
}


//!Left Shift with shift off
//!Shifts the value to the left by n bits.
//! The bits that get shifted off val are returned in lsb of return value
//!\param [in,out] val  Value to shift
//!\param [in] old_val Values to shift on to val
//!\returns value with bits that are shifted off val in lsb
template <typename T>
T lsoso(T& val, T& old_val, uint8_t n) {
	T temp(val);
	val = (val << n) | old_val;
    return temp >> ((sizeof(T)*CHAR_BIT) - n);
};

//!Right Shift with shift off
//!Shifts the value to the right by n bits.
//! The bits that would get shifted off value are returned in msb of return value
//!\param [in,out] val  Value to shift
//!\param [in] old_val Values to shift on to val
template <typename T>
T rsoso(T& val, T& old_val, uint8_t n) {
	T temp(val);
	val = (val>>n) | old_val;
    return temp << ((sizeof(T)*CHAR_BIT) - n);
}



#endif /* defined(__DynamicBitset__bitwise_ops__) */
