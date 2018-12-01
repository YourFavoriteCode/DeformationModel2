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
	

	/********************************************************
	********		Общие параметры программы		*********
	********************************************************/

	// Отладочный режим
	bool CONFIG_DEBUG						= false;
	// Кол-во вычислительных потоков OMP
	int ompThreadCount						= 4;
	

	/********************************************************
	**********		   Прочие параметры			 ************
	********************************************************/

	// Симметризация диад nb
	bool isSymmetrycal						= true;
	// Считывание ориентаций и нормалей из файлы
	int fixedOrientations					= 0;
	// Шаг интегрирования
	double dt								= 5e-3;
	// Считывание остаточных напряжений из файла и запись в конце
	bool usingInititalStress				= false;
	// Сохранение информации о ССТ
	bool usingStandardTriangleSaving		= false;
	bool usingFragmentation					= true;
	double fragmentationCriteria			= 4.5e7;
	

	/********************************************************
	**********		Параметры поликристалла	     ************
	********************************************************/

	// Используемый материал
	int materialType = 0;
	// Кол-во зерен
	int grainCountLinear = 2;
	// Закон распределения размеров зерен
	DistributionType grainSizeDistribLaw = DISTRIB_NORMAL;
	// Мат.ожидание размера зерна
	double grainSizeDistribM = 5e-5;
	// Дисперсия размера зерна
	double grainSizeDistribD = 0;
	// Кол-во учитываемых соседей зерна DEPRECATED
	int grainSurroundCount = 6;
	// Процент основной фазы в материале
	int mainPhasePercent = 100;
	// Случайное распределение начальных ориентаций
	bool randomOrientations = true;
	// Тип задания ориентации
	int orientationType = 0;


	/********************************************************
	**********   Упруговязкопластический закон   ************
	********************************************************/

	// Начальная скорость сдвига при критическом касательном напряжении
	double shearRateLawDgm0 = 1e-5;
	// Степенной параметр скоростной чувствительности
	double shearRateLawM = 83;


	/********************************************************
	**********   Различные параметры нагружения   ***********
	********************************************************/

	// Условие одноосности НДС
	bool trueUniaxial = false;
	// Упругая разгрузка после каждого цикла нагружения
	bool withUnloading = false;
	// Кол-во циклов нагружения
	int loadCycleCount = 1;
	// Предел интенсивности деформаций
	double maxStrainIntencity = 2e-1;
	// Градиент места
	model::Tensor gradV(0, 0, -0.0015, 0, -0.0015, 0, 0, 0, 0.003);


	/********************************************************
	*****************       Ротации      ********************
	********************************************************/

	// Модель ротаций Тейлора
	bool usingRotationsTaylor				= false;
	// Модель несовместности свдигов
	bool usingRotationsTrusov				= false;
	// Ротационное упрочнение
	bool usingRotationsHardening			= false;

	// Параметр упругой составляющей разворота
	double rotationParamA					= 3e-9;
	// Параметр пластической составляющей
	double rotationParamH					= 1e-7;
	// Параметр лямбда
	double rotationParamL					= 10;
	// Начальный критический момент
	double rotationParamMc					= 3e4;
	// Параметры ротационного упрочнения
	double rotationParamHardK1;
	double rotationParamHardK2;
	

	/********************************************************
	*****************     Упрочнения     ********************
	********************************************************/
	
	// Базовое слагаемое упрочнения
	bool usingHardeningBase					= false;
	// Зерно-граничное упрочнение
	bool usingHardeningBound				= false;

	// Параметр зернограничного упрочнения
	double hardeningParamBoundK				= 5e6;
	// Базовое слагаемое: парметр дельта
	double hardeningParamBaseDelta			= -0.2;
	// Базовое слагаемое: параметр пси
	double hardeningParamBasePsi			= 1.00002;
	// Базовое слагаемое: парметр А
	double hardeningParamBaseA				= 0.0057;


	/********************************************************
	***********			Выходные файлы		*****************
	********************************************************/

	// Сохранение интенсивностей тензоров
	bool saveIntensity						= true;
	// Сохранение данных макроуровня
	bool saveMacro							= true;
	// Сохранение данных мехоуровня
	bool saveMeso							= false;
	// Сохранение активных СС
	bool saveActiveSS						= false;

	// Период сохранения диаграммы НДС (в % выполнения)
	double periodSavePlot = 2;
	// Период сохранения ПФ (в % выполнения)
	double periodSavePolus = 25;
	// Период сохранения отладочных данных
	int saveVariablesPeriodStep = 0;
	// Нижнее ограничение записи отладочной информации
	int saveVariablesStartStep = 0;
	// Верхнее ограничение записи отладочной информации
	int saveVariablesStopStep = INT_MAX;

	// Покомпонентное сохранение тензоров
	bool save11								= false;
	bool save12								= false;
	bool save13								= false;
	bool save21								= false;
	bool save22								= false;
	bool save23								= false;
	bool save31								= false;
	bool save32								= false;
	bool save33								= false;


	// Чтение целочисленных параметров
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