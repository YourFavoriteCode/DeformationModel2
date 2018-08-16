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
	bool isSymmetrycal					= true;
	bool trueUniaxial					= false;
	bool withUnloading					= false;
	bool randomOrientations				= true;
	int orientationType					= 0;
	int fixedOrientations				= 0;
	double dt							= 5e-4;
	int materialType					= 0;
	double maxStrainIntencity			= 1e-1;
	model::Tensor gradV;
	
	int surroundCount					= 6;
	int material_purity					= 100;
	bool read_init_stress				= false;
	
	int fragm_size_law					= 0;
	double fragm_size_m					= 5e-5;
	double fragm_size_dsp				= 0;

	int fragm_count						= 4;
	int cycle_count						= 1;
	int thread_count					= 1;

	double plot_period					= 2;
	double polus_period					= 25;
	int debug_period					= 0;
	int DEBUG_START						= 0;
	int DEBUG_STOP						= INT_MAX;

	double dgm0							= 1e-5;
	double m							= 100;

	bool ROTATIONS_TAYLOR				= false;
	bool ROTATIONS_TRUSOV				= false;
	bool ROTATIONS_HARDENING			= false;

	double rotationParamA				= 3e-9;
	double rotationParamH				= 1e-7;
	double rotationParamL				= 10;
	double rotationParamMc				= 3e4;

	double rotationParamHardK1;
	double rotationParamHardK2;

	bool usingHardeningBase				= false;
	bool usingHardeningBound			= false;

	double hardeningParamBoundK			= 0;
	double hardeningParamBaseDelta		= 0;
	double hardeningParamBasePsi		= 0;
	double hardeningParamBaseA			= 0;

	int SurroundsGrade					= 1;
	bool SST_SAVING						= false;

	bool FRAGMENTATION					= false;
	int Grain_size						= 3;

	
	bool isSaveIntensity				= true;
	bool isSaveMacro					= true;
	bool isSaveMeso						= false;
	bool isSaveActiveSS					= false;

	bool save11							= false;
	bool save12							= false;
	bool save13							= false;
	bool save21							= false;
	bool save22							= false;
	bool save23							= false;
	bool save31							= false;
	bool save32							= false;
	bool save33							= false;

	int getValue(tinyxml2::XMLElement *rootnode, const char* name, int* var)//Целочисленные
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
	
	int getValue(tinyxml2::XMLElement *rootnode, const char* name, double* var)//Вещественные
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

	int getValue(tinyxml2::XMLElement *rootnode, const char* name, bool* var)//Логические
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
		getValue(rootnode, "GrainCount", &fragm_count);
		getValue(rootnode, "Material", &materialType);
		getValue(rootnode, "RandomOrientations", &randomOrientations);
		getValue(rootnode, "M", &m);
		getValue(rootnode, "dGamma0", &dgm0);
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
		int rt;
		getValue(rootnode, "Rotate", &rt);
		if (rt == 0) ROTATIONS_TAYLOR = true;
		else if (rt == 1) ROTATIONS_TRUSOV = true;
		getValue(rootnode, "RotateA", &rotationParamA);
		getValue(rootnode, "RotateH", &rotationParamH);
		getValue(rootnode, "RotateLambda", &rotationParamL);
		getValue(rootnode, "RotateMc", &rotationParamMc);

		getValue(rootnode, "CycleCount", &cycle_count);
		getValue(rootnode, "Unloading", &withUnloading);
		getValue(rootnode, "SurroundsDegree", &SurroundsGrade);
		getValue(rootnode, "RotationHardening", &ROTATIONS_HARDENING);
		getValue(rootnode, "MaterialPurity", &material_purity);
		getValue(rootnode, "PlotPeriod", &plot_period);
		getValue(rootnode, "PolusPeriod", &polus_period);
		getValue(rootnode, "DebugPeriod", &debug_period);
		getValue(rootnode, "ThreadCount", &thread_count);
		getValue(rootnode, "FixedOrientations", &fixedOrientations);
		getValue(rootnode, "OrientationType", &orientationType);

		/*Градиент скорости*/
		gradV.setZero();
		getValue(rootnode, "gradV00", &gradV.C[0][0]);
		getValue(rootnode, "gradV01", &gradV.C[0][1]);
		getValue(rootnode, "gradV02", &gradV.C[0][2]);
		getValue(rootnode, "gradV10", &gradV.C[1][0]);
		getValue(rootnode, "gradV11", &gradV.C[1][1]);
		getValue(rootnode, "gradV12", &gradV.C[1][2]);
		getValue(rootnode, "gradV20", &gradV.C[2][0]);
		getValue(rootnode, "gradV21", &gradV.C[2][1]);
		getValue(rootnode, "gradV22", &gradV.C[2][2]);
	
		getValue(rootnode, "FragmSizeLaw", &fragm_size_law);
		getValue(rootnode, "FragmSizeM", &fragm_size_m);
		getValue(rootnode, "FragmSizeDsp", &fragm_size_dsp);
		getValue(rootnode, "StartWritingDbgInfo", &DEBUG_START);
		getValue(rootnode, "StopWritingDbgInfo", &DEBUG_STOP);

		/*Выходные данные*/
		getValue(rootnode, "SaveIntense", &isSaveIntensity);
		getValue(rootnode, "SaveMacro", &isSaveMacro);
		getValue(rootnode, "SaveMeso", &isSaveMeso);
		getValue(rootnode, "SaveActiveSS", &isSaveActiveSS);

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
		getValue(rootnode, "ReadInitStress", &read_init_stress);
		getValue(rootnode, "SaveSST", &SST_SAVING);
		getValue(rootnode, "Fragmentation", &FRAGMENTATION);
		getValue(rootnode, "GrainSize", &Grain_size);

		getValue(rootnode, "ROT_HARD_K1", &rotationParamHardK1);
		getValue(rootnode, "ROT_HARD_K2", &rotationParamHardK2);
	

		return 0;
	}
}