#pragma once
#ifndef __PARAMS_H 
#define __PARAMS_H

#include "Tensor.h"

namespace prms
{
	/********************************************************
	**********   Различные параметры нагружения   ***********
	********************************************************/
	extern bool trueUniaxial;					//Условие одноосности
	extern bool withUnloading;					//Упругая разгрузка
	extern double maxStrainIntencity;				//Предел интенсивности деформаций
	extern int cycle_count;					//Кол-во циклов нагружения
	extern model::Tensor gradV;					//Градиент места

	/********************************************************
	**********		   Прочие параметры			 ************
	********************************************************/
	extern bool isSymmetrycal;					//Симметризация диад nb
	extern bool randomOrientations;				//Случайное распределение начальных ориентаций
	extern int orientationType;					//Тип задания ориентации
	extern int fixedOrientations;					//Считывание ориентаций и нормалей из файлы
	extern double dt;						//Шаг интегрирования
	extern int thread_count;				//Кол-во потоков
	extern int debug_period;				//Период сохранения отладочных данных
	extern int DEBUG_START;					//Нижнее ограничение записи отладочной информации
	extern int DEBUG_STOP;					//Верхнее ограничение записи отладочной информации
	extern bool read_init_stress;			//Считывание остаточных напряжений из файла и запись в конце
	extern bool SST_SAVING;					//Сохранение информации о ССТ

	/********************************************************
	**********		Параметры поликристалла	     ************
	********************************************************/
	extern int materialType;					//Используемый материал
	extern int SurroundsGrade;				//Степень учёта соседних элементов
	extern int surroundCount;				//Кол-во учитываемых соседей
	extern int fragm_count;					//Кол-во фрагментов
	extern int material_purity;				//Процент основной фазы в материале
	
	extern int fragm_size_law;				//Закон распределения размеров фрагментов
	extern double fragm_size_m;				//Мат.ожидание размера зерна
	extern double fragm_size_dsp;			//Дисперсия размера зерна

	/********************************************************
	**********   Упруговязкопластический закон   ************
	********************************************************/

	extern double dgm0;						//Начальная скорость сдвига при критическом касательном напряжении
	extern double m;						//Степенной параметр скоростной чувствительности

	/********************************************************
	*****************       Ротации      ********************
	********************************************************/

	extern bool ROTATIONS_TAYLOR;			//Модель ротаций Тейлора
	extern bool ROTATIONS_TRUSOV;			//Модель несовместности свдигов
	extern bool ROTATIONS_HARDENING;		//Ротационное упрочнение

	extern double rotationParamA;					//Параметр упругой составляющей
	extern double rotationParamH;					//Параметр пластической составляющей
	extern double rotationParamL;					//Параметр лямбда
	extern double rotationParamMc;					//Начальный критический момент

	extern double rotationParamHardK1;
	extern double rotationParamHardK2;
	
	/********************************************************
	*****************     Упрочнения     ********************
	********************************************************/

	extern bool usingHardeningBase;				//Базовое слагаемое упрочнения
	extern bool usingHardeningBound;			//Зерно-граничное упрочнение

	extern double hardeningParamBoundK;				//Параметр зернограничного упрочнения
	extern double hardeningParamBaseDelta;			//Парметр дельта
	extern double hardeningParamBasePsi;			//Параметр пси
	extern double hardeningParamBaseA;				//Парметр А
	
	/********************************************************
	***********			Фрагментация		*****************
	********************************************************/

	extern bool FRAGMENTATION;
	extern int Grain_size;

	/********************************************************
	***********			Выходные файлы		*****************
	********************************************************/

	extern double plot_period;				//Период сохранения диаграммы НДС
	extern double polus_period;				//Период сохранения ПФ
	extern bool isSaveIntensity;				//Сохранение интенсивностей тензоров
	extern bool isSaveMacro;					//Сохранение данных макроуровня
	extern bool isSaveMeso;					//Сохранение данных мехоуровня
	extern bool isSaveActiveSS;				//Сохранение активных СС	
	extern bool save11;						//Покомпонентное сохранение
	extern bool save12;
	extern bool save13;
	extern bool save21;
	extern bool save22;
	extern bool save23;
	extern bool save31;
	extern bool save32;
	extern bool save33;

	int ReadParams(const char *);			//Считывание параметров из файла

}

#endif __PARAMS_H
