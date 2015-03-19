/*
 *  statistics.h
 *
 *  Created by Paul Lott on 4/22/09.
 *  Copyright 2009 University of California, Davis. All rights reserved.
 *
 */

#ifndef DISTRIBUTIONS_H_
#define DISTRIBUTIONS_H_

#include <vector>
#include <math.h>
#include "stochMath.h"
#include <iostream>


namespace StochHMM {
    
#define PI 3.141592653589793238463
    
        
    ///////////////////DISTRIBUTIONS////////////////////////
    
    
    //Discrete & Finite
    void sBinomial (std::vector<double>&, int, double);
    void cBinomial (std::vector<double>&, int, double);
    
    
    void sBetaBinomial (std::vector<double>&, int, double, double);
    void cBetaBinomial (std::vector<double>&, int, double, double);
    
    
    void sDegenerate(std::vector<double>&, double);
    void cDegenerate(std::vector<double>&, double);
    
    void sUniform (std::vector<double> &,int,int);
    void cUniform (std::vector<double> &,int,int);
    
    void sHypergeometric (std::vector<double>&,int, int, int);
    void cHypergeometric (std::vector<double>&,int, int, int);
    
    //void sPoissonBinomial (std::vector<double>&, int, double);
    //void cPoissonBinomial (std::vector<double>&, int, double);
    
    //Continuous
    
    void sCauchy(std::vector<double>&,double, double);
    
    
    void sChiSquared(std::vector<double>&, double);
    double cChiSquared(double x, double df);
    
    
    void sExponential (std::vector<double>&, double);
    void sExtremeValue(std::vector<double>&, double, double,double);
    void sFDistribution(std::vector<double>&, double, double);
    void sGamma(std::vector<double>&, double, double);
    void sGeometric (std::vector<double>&, double);
    void sLaplace(std::vector<double>&, double, double);
    void sLogNormal(std::vector<double>&,double, double);
    void sLogarithmic(std::vector<double>&, double);
    void sLogistic(std::vector<double>&, double, double);
    void sMaxwellBoltzman(std::vector<double>&, double);
    void sNegativeBinomial(std::vector<double>&,int, double);
    void sNonCentralF(std::vector<double>&,double, double, double);
    void sNonCentralT(std::vector<double>&, double, double);
    void sNormal(std::vector<double>&, double, double);
    void sPareto(std::vector<double>&,double, double);
    void sPoisson(std::vector<double>&, double);
    void sRayleigh(std::vector<double>&, double);
    void sTriangular(std::vector<double>&,int, int, int);
    void sUser (std::vector<double>&, std::vector<double>&);
    void sWeibull(std::vector<double>&,double, double);
}



#endif /*INC_FILENAME_H*/
