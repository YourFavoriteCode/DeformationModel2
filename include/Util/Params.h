#pragma once
#ifndef __PARAMS_H 
#define __PARAMS_H

#include "Tensor.h"

namespace prms
{

	// Типы распределения
	enum DistributionType {
		DISTRIB_UNIFORM,	// Равномерное
		DISTRIB_NORMAL,		// Нормальное
		DISTRIB_LOGNORMAL,	// Логнормальное
		DISTRIB_EXPONENT	// Показательное
	};


	
	
	extern bool CONFIG_DEBUG;					


	extern bool trueUniaxial;
	extern bool withUnloading;
	extern double maxStrainIntencity;
	extern int loadCycleCount;
	extern model::Tensor gradV;


	extern bool isSymmetrycal;
	extern bool randomOrientations;
	extern int orientationType;
	extern int fixedOrientations;
	extern double dt;
	extern int ompThreadCount;
	extern int saveVariablesPeriodStep;
	extern int saveVariablesStartStep;
	extern int saveVariablesStopStep;
	extern bool usingInititalStress;
	extern bool usingStandardTriangleSaving;


	extern int materialType;
	extern int grainSurroundCount;
	extern int grainCountLinear;
	extern int mainPhasePercent;
	extern DistributionType grainSizeDistribLaw;
	extern double grainSizeDistribM;
	extern double grainSizeDistribD;

	extern double shearRateLawDgm0;
	extern double shearRateLawM;

	extern bool usingRotationsTaylor;			
	extern bool usingRotationsTrusov;			
	extern bool usingRotationsHardening;		

	extern double rotationParamA;
	extern double rotationParamH;
	extern double rotationParamL;
	extern double rotationParamMc;

	extern double rotationParamHardK1;
	extern double rotationParamHardK2;
	
	extern bool usingHardeningBase;
	extern bool usingHardeningBound;

	extern double hardeningParamBoundK;
	extern double hardeningParamBaseDelta;
	extern double hardeningParamBasePsi;
	extern double hardeningParamBaseA;

	extern double periodSavePlot;
	extern double periodSavePolus;
	extern bool saveIntensity;
	extern bool saveMacro;
	extern bool saveMeso;
	extern bool saveActiveSS;
	
	extern bool save11;						
	extern bool save12;
	extern bool save13;
	extern bool save21;
	extern bool save22;
	extern bool save23;
	extern bool save31;
	extern bool save32;
	extern bool save33;

	extern bool usingFragmentation;
	extern double fragmentationCriteria;

	//Считывание параметров из файла
	int ReadParams(const char *);			

}

#endif __PARAMS_H
