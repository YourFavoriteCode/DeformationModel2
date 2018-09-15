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

	/********************************************************
	**********   Различные параметры нагружения   ***********
	********************************************************/
	extern bool trueUniaxial;					//Условие одноосности
	extern bool withUnloading;					//Упругая разгрузка
	extern double maxStrainIntencity;			//Предел интенсивности деформаций
	extern int loadCycleCount;					//Кол-во циклов нагружения
	extern model::Tensor gradV;					//Градиент места

	/********************************************************
	**********		   Прочие параметры			 ************
	********************************************************/
	extern bool isSymmetrycal;					//Симметризация диад nb
	extern bool randomOrientations;				//Случайное распределение начальных ориентаций
	extern int orientationType;					//Тип задания ориентации
	extern int fixedOrientations;				//Считывание ориентаций и нормалей из файлы
	extern double dt;							//Шаг интегрирования
	extern int ompThreadCount;					//Кол-во потоков
	extern int saveVariablesPeriodStep;			//Период сохранения отладочных данных
	extern int saveVariablesStartStep;			//Нижнее ограничение записи отладочной информации
	extern int saveVariablesStopStep;			//Верхнее ограничение записи отладочной информации
	extern bool usingInititalStress;			//Считывание остаточных напряжений из файла и запись в конце
	extern bool usingStandardTriangleSaving;	//Сохранение информации о ССТ

	/********************************************************
	**********		Параметры поликристалла	     ************
	********************************************************/
	extern int materialType;					//Используемый материал
	extern int grainSurroundCount;				//Кол-во учитываемых соседей
	extern int grainCountLinear;				//Кол-во фрагментов
	extern int mainPhasePercent;				//Процент основной фазы в материале
	extern DistributionType grainSizeDistribLaw;//Закон распределения размеров фрагментов
	extern double grainSizeDistribM;			//Мат.ожидание размера зерна
	extern double grainSizeDistribD;			//Дисперсия размера зерна

	/********************************************************
	**********   Упруговязкопластический закон   ************
	********************************************************/

	extern double shearRateLawDgm0;				//Начальная скорость сдвига при критическом касательном напряжении
	extern double shearRateLawM;				//Степенной параметр скоростной чувствительности

	/********************************************************
	*****************       Ротации      ********************
	********************************************************/

	extern bool usingRotationsTaylor;			//Модель ротаций Тейлора
	extern bool usingRotationsTrusov;			//Модель несовместности свдигов
	extern bool usingRotationsHardening;		//Ротационное упрочнение

	extern double rotationParamA;				//Параметр упругой составляющей
	extern double rotationParamH;				//Параметр пластической составляющей
	extern double rotationParamL;				//Параметр лямбда
	extern double rotationParamMc;				//Начальный критический момент

	extern double rotationParamHardK1;			//Параметры ротационного упрочнения
	extern double rotationParamHardK2;
	
	/********************************************************
	*****************     Упрочнения     ********************
	********************************************************/

	extern bool usingHardeningBase;				//Базовое слагаемое упрочнения
	extern bool usingHardeningBound;			//Зерно-граничное упрочнение

	extern double hardeningParamBoundK;			//Параметр зернограничного упрочнения
	extern double hardeningParamBaseDelta;		//Парметр дельта
	extern double hardeningParamBasePsi;		//Параметр пси
	extern double hardeningParamBaseA;			//Парметр А

	/********************************************************
	***********			Выходные файлы		*****************
	********************************************************/

	extern double periodSavePlot;			//Период сохранения диаграммы НДС
	extern double periodSavePolus;			//Период сохранения ПФ
	extern bool saveIntensity;				//Сохранение интенсивностей тензоров
	extern bool saveMacro;					//Сохранение данных макроуровня
	extern bool saveMeso;					//Сохранение данных мехоуровня
	extern bool saveActiveSS;				//Сохранение активных СС	
	
	//Покомпонентное сохранение тензоров

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

	//Считывание параметров из файла
	int ReadParams(const char *);			

}

#endif __PARAMS_H
