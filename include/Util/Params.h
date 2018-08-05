﻿#pragma once
#ifndef __PARAMS_H 
#define __PARAMS_H

#include "Tensor.h"

namespace prms
{
	/********************************************************
	**********   Различные параметры нагружения   ***********
	********************************************************/
	extern bool REAL_UNIAX;					//Условие одноосности
	extern bool UNLOADING;					//Упругая разгрузка
	extern double strain_max;				//Предел интенсивности деформаций
	extern int cycle_count;					//Кол-во циклов нагружения
	extern model::Tensor gradV;					//Градиент места

	/********************************************************
	**********		   Прочие параметры			 ************
	********************************************************/
	extern bool SYMMETRY;					//Симметризация диад nb
	extern bool RAND_ORIENT;				//Случайное распределение начальных ориентаций
	extern int ORIENT_TYPE;					//Тип задания ориентации
	extern int fix_orient;					//Считывание ориентаций и нормалей из файлы
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
	extern int material;					//Используемый материал
	extern int SurroundsGrade;				//Степень учёта соседних элементов
	extern int surround_count;				//Кол-во учитываемых соседей
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

	extern double ROT_A;					//Параметр упругой составляющей
	extern double ROT_H;					//Параметр пластической составляющей
	extern double ROT_L;					//Параметр лямбда
	extern double ROT_MC;					//Начальный критический момент

	extern double ROT_HARD_K1;
	extern double ROT_HARD_K2;
	
	/********************************************************
	*****************     Упрочнения     ********************
	********************************************************/

	extern bool HARDENING_BASE;				//Базовое слагаемое упрочнения
	extern bool HARDENING_BOUND;			//Зерно-граничное упрочнение

	extern double HARD_BOUND_K;				//Параметр зернограничного упрочнения
	extern double HARD_BASE_DELTA;			//Парметр дельта
	extern double HARD_BASE_PSI;			//Параметр пси
	extern double HARD_BASE_A;				//Парметр А
	
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
	extern bool SaveIntense;				//Сохранение интенсивностей тензоров
	extern bool SaveMacro;					//Сохранение данных макроуровня
	extern bool SaveMeso;					//Сохранение данных мехоуровня
	extern bool SaveActiveSS;				//Сохранение активных СС	
	extern bool Save11;						//Покомпонентное сохранение
	extern bool Save12;
	extern bool Save13;
	extern bool Save21;
	extern bool Save22;
	extern bool Save23;
	extern bool Save31;
	extern bool Save32;
	extern bool Save33;

	int ReadParams(const char *);			//Считывание параметров из файла

}

#endif __PARAMS_H
