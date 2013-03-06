//
//  bitwise_ops.cpp
//  DynamicBitset
//
//  Created by Paul Lott on 3/1/13.
//  Copyright (c) 2013 Korf Lab, Genome Center, UC Davis, Davis, CA. All rights reserved.
//

#include "bitwise_ops.h"

std::string intToBinString(uint32_t val){
	std::string output("");
	for(size_t i = 31; i != SIZE_MAX ;--i){
		if (val & (1 << i)){
			output+="1";
		}
		else{
			output+="0";
		}
	}
	return output;
}

//Finds most significant bit set
size_t msb(uint32_t val){
	const uint32_t mask[] = {
		0x00000FFFF,
		0x0000000FF,
		0x00000000F,
		0x000000003,
		0x000000001
	};
	int hi = 32;
	int lo = 0;
	int i = 0;
	
	if (val == 0)
		return 0;
	
	for (i = 0; i < 5; i++) {
		int mi = lo + (hi - lo) / 2;
		
		if ((val >> mi) != 0)
			lo = mi;
		else if ((val & (mask[i] << lo)) != 0)
			hi = mi;
	}
	
	return lo;
}


//Find lowest lsb 
size_t lowestBitSet(uint32_t x){
	static const unsigned int MultiplyDeBruijnBitPosition[32] =
	{
		// precomputed lookup table
		0,  1, 28,  2, 29, 14, 24,  3, 30, 22, 20, 15, 25, 17,  4,  8,
		31, 27, 13, 23, 21, 19, 16,  7, 26, 12, 18,  6, 11,  5, 10,  9
	};
	
	// leave only lowest bit
	x  &= -int(x);
	// DeBruijn constant
	x  *= 0x077CB531;
	// get upper 5 bits
	x >>= 27;
	// convert to actual position
	return MultiplyDeBruijnBitPosition[x];
}

//!Population distance (Hamming Distance)
//!Returns the number of bits set in the 32bit integer
int8_t popCount(uint32_t val){
	val = val - ((val >> 1) & 0x55555555);                    // reuse input as temporary
	val = (val & 0x33333333) + ((val >> 2) & 0x33333333);     // temp
	return ((val + (val >> 4) & 0xF0F0F0F) * 0x1010101) >> 24; // count
}



//!Returns integer's higher bits from position up to msb
//!Lower bits are set to zero
uint32_t maskLeft(uint32_t val,uint8_t pos){
	if (pos>=32){
		return 0;
	}
	
	return val & ~((1 << pos) - 1);
}

//!Returns integer's lower bits from position down to lsb
//!Higher bits are set to zero
uint32_t maskRight(uint32_t val, uint8_t pos){
	if (pos >31){
		return 0;
	}
	else if (pos == 31){
		return val & (1<<31);
	}	
	return val & ((1 << pos+1) - 1);
}

////Rotate left with carry through
//uint32_t ror(uint32_t& val, uint8_t shift_register) {
//	uint32_t temp(val);
//	//Shifted integer
//	val = val << shift_register;
//	
//	//Bits shifted off register
//	return temp >> (32 - shift_register);
//    //return (val << shift_register) | (val >> (32 - shift_register));
//}
//
////Rotate right with carry through
//uint32_t ror(uint32_t& val, uint8_t shift_register) {
//	uint32_t temp(val);
//	//Shifted integer
//	val = val >> shift_register;
//	
//	//Bits shifted off register
//	return temp >> (32 - shift_register);
//    //return (val >> shift_register) | (val << (32 - shift_register));
//}






//unsigned char rol (unsigned char x)
//{
//	return ((x >> 7) | (x << 1));
//}
//
//unsigned char ror (unsigned char x)
//{
//	return ((x << 7) | (x >> 1));
//}
