// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"

#include "Params.h"
#include "tinyxml2.h"

/*
Параметры модели, доступные из любого места программы
*/

namespace prms
{
	/*****************************************************************
	*****		Задание значений параметров по-умолчанию		******
	*****														******
	*****	Если в файле параметров не будет найден параметр,	******
	*****	его значение будет взято по-умолчанию				******
	*****************************************************************/
	bool CONFIG_DEBUG						= false;
	
	bool isSymmetrycal						= true;
	bool trueUniaxial						= false;
	bool withUnloading						= false;
	bool randomOrientations					= true;
	int orientationType						= 0;
	int fixedOrientations					= 0;
	double dt								= 5e-3;
	int materialType						= 0;
	double maxStrainIntencity				= 1e-1;
	model::Tensor gradV;
	
	int grainSurroundCount					= 6;
	int mainPhasePercent					= 100;
	bool usingInititalStress				= false;

	DistributionType grainSizeDistribLaw	= DISTRIB_NORMAL;
	double grainSizeDistribM				= 5e-5;
	double grainSizeDistribD				= 0;

	int grainCountLinear					= 2;
	int loadCycleCount						= 1;
	int ompThreadCount						= 1;

	double periodSavePlot					= 2;
	double periodSavePolus					= 25;
	int saveVariablesPeriodStep				= 0;
	int saveVariablesStartStep				= 0;
	int saveVariablesStopStep				= INT_MAX;

	double shearRateLawDgm0					= 1e-5;
	double shearRateLawM					= 100;

	bool usingRotationsTaylor				= false;
	bool usingRotationsTrusov				= false;
	bool usingRotationsHardening			= false;

	double rotationParamA					= 3e-9;
	double rotationParamH					= 1e-7;
	double rotationParamL					= 10;
	double rotationParamMc					= 3e4;

	double rotationParamHardK1;
	double rotationParamHardK2;

	bool usingHardeningBase					= false;
	bool usingHardeningBound				= false;

	double hardeningParamBoundK				= 0;
	double hardeningParamBaseDelta			= 0;
	double hardeningParamBasePsi			= 0;
	double hardeningParamBaseA				= 0;

	bool usingFragmentation					= true;
	double fragmentationCriteria			= 1e10;

	bool usingStandardTriangleSaving		= false;

	bool saveIntensity						= true;
	bool saveMacro							= true;
	bool saveMeso							= false;
	bool saveActiveSS						= false;

	bool save11								= false;
	bool save12								= false;
	bool save13								= false;
	bool save21								= false;
	bool save22								= false;
	bool save23								= false;
	bool save31								= false;
	bool save32								= false;
	bool save33								= false;

	//Чтение целочисленных параметров
	int getValue(tinyxml2::XMLElement *rootnode, const char* name, int* var)
	{
		if (rootnode->FirstChildElement(name) != NULL)//Проверка, есть ли в файле такой элемент
		{
			const char* title;
			title = rootnode->FirstChildElement(name)->GetText();
			*var = atoi(title);
		}
		else//Если нет, то предупреждаем и выводим значение по умолчанию
		{
			printf(" Parameter <%s> not found. Default value: %d.\n", name, *var);
		}
		return 1;
	}
	
	//Чтение вещественных параметров
	int getValue(tinyxml2::XMLElement *rootnode, const char* name, double* var)
	{
		if (rootnode->FirstChildElement(name) != NULL)//Проверка, есть ли в файле такой элемент
		{
			const char* title;
			title = rootnode->FirstChildElement(name)->GetText();
			*var = atof(title);
		}
		else//Если нет, то предупреждаем и выводим значение по умолчанию
		{
			printf(" Parameter <%s> not found. Default value: %e.\n", name, *var);
		}
		return 1;
	}

	//Чтение логических параметров
	int getValue(tinyxml2::XMLElement *rootnode, const char* name, bool* var)
	{
		if (rootnode->FirstChildElement(name) != NULL)//Проверка, есть ли в файле такой элемент
		{
			const char* title;
			title = rootnode->FirstChildElement(name)->GetText();
			*var = (bool)atoi(title);
		}
		else//Если нет, то предупреждаем и выводим значение по умолчанию
		{
			printf(" Parameter <%s> not found. Default value: %d.\n", name, *var);
		}
		return 1;
	}

	int ReadParams(const char * filename)
	{
		tinyxml2::XMLDocument doc;
		doc.LoadFile(filename);//Загрузка файла

		tinyxml2::XMLElement *rootnode = doc.FirstChildElement("Parameters");
		if (rootnode == NULL) return 1;

		getValue(rootnode, "GrainCount", &grainCountLinear);
		getValue(rootnode, "Material", &materialType);
		getValue(rootnode, "RandomOrientations", &randomOrientations);
		getValue(rootnode, "M", &shearRateLawM);
		getValue(rootnode, "dGamma0", &shearRateLawDgm0);
		getValue(rootnode, "DeformationLimit", &maxStrainIntencity);
		getValue(rootnode, "IntegrationStep", &dt);
		getValue(rootnode, "RealUniaxial", &trueUniaxial);
		getValue(rootnode, "Symmetrisation", &isSymmetrycal);
		/*Упрочнение*/
		getValue(rootnode, "HardBase", &usingHardeningBase);
		getValue(rootnode, "HardBound", &usingHardeningBound);
		getValue(rootnode, "HardBaseDelta", &hardeningParamBaseDelta);
		getValue(rootnode, "HardBasePsi", &hardeningParamBasePsi);
		getValue(rootnode, "HardBaseA", &hardeningParamBaseA);
		getValue(rootnode, "HardBoundK", &hardeningParamBoundK);
		/*Ротации*/
		int buf;
		getValue(rootnode, "Rotate", &buf);
		if (buf == 0) usingRotationsTaylor = true;
		else if (buf == 1) usingRotationsTrusov = true;
		getValue(rootnode, "RotateA", &rotationParamA);
		getValue(rootnode, "RotateH", &rotationParamH);
		getValue(rootnode, "RotateLambda", &rotationParamL);
		getValue(rootnode, "RotateMc", &rotationParamMc);

		getValue(rootnode, "CycleCount", &loadCycleCount);
		getValue(rootnode, "Unloading", &withUnloading);
		getValue(rootnode, "RotationHardening", &usingRotationsHardening);
		getValue(rootnode, "MaterialPurity", &mainPhasePercent);
		getValue(rootnode, "PlotPeriod", &periodSavePlot);
		getValue(rootnode, "PolusPeriod", &periodSavePolus);
		getValue(rootnode, "DebugPeriod", &saveVariablesPeriodStep);
		getValue(rootnode, "ThreadCount", &ompThreadCount);
		getValue(rootnode, "FixedOrientations", &fixedOrientations);
		getValue(rootnode, "OrientationType", &orientationType);

		/*Градиент скорости*/
		gradV.setZero();
		getValue(rootnode, "gradV00", &gradV.c[0][0]);
		getValue(rootnode, "gradV01", &gradV.c[0][1]);
		getValue(rootnode, "gradV02", &gradV.c[0][2]);
		getValue(rootnode, "gradV10", &gradV.c[1][0]);
		getValue(rootnode, "gradV11", &gradV.c[1][1]);
		getValue(rootnode, "gradV12", &gradV.c[1][2]);
		getValue(rootnode, "gradV20", &gradV.c[2][0]);
		getValue(rootnode, "gradV21", &gradV.c[2][1]);
		getValue(rootnode, "gradV22", &gradV.c[2][2]);
	
		getValue(rootnode, "FragmSizeLaw", &buf);
		grainSizeDistribLaw = DistributionType(buf);
		getValue(rootnode, "FragmSizeM", &grainSizeDistribM);
		getValue(rootnode, "FragmSizeDsp", &grainSizeDistribD);
		getValue(rootnode, "StartWritingDbgInfo", &saveVariablesStartStep);
		getValue(rootnode, "StopWritingDbgInfo", &saveVariablesStopStep);

		/*Выходные данные*/
		getValue(rootnode, "SaveIntense", &saveIntensity);
		getValue(rootnode, "SaveMacro", &saveMacro);
		getValue(rootnode, "SaveMeso", &saveMeso);
		getValue(rootnode, "SaveActiveSS", &saveActiveSS);

		getValue(rootnode, "Save11", &save11);
		getValue(rootnode, "Save12", &save12);
		getValue(rootnode, "Save13", &save13);
		getValue(rootnode, "Save21", &save21);
		getValue(rootnode, "Save22", &save22);
		getValue(rootnode, "Save23", &save23);
		getValue(rootnode, "Save31", &save31);
		getValue(rootnode, "Save32", &save32);
		getValue(rootnode, "Save33", &save33);

		/*Экспериментальные параметры*/
		getValue(rootnode, "ReadInitStress", &usingInititalStress);
		getValue(rootnode, "SaveSST", &usingStandardTriangleSaving);

		getValue(rootnode, "ROT_HARD_K1", &rotationParamHardK1);
		getValue(rootnode, "ROT_HARD_K2", &rotationParamHardK2);

		getValue(rootnode, "using_fragmentation", &usingFragmentation);
		getValue(rootnode, "fragmentation_criteria", &fragmentationCriteria);
		return 0;
	}
}