#pragma once
#ifndef __POLYCRYST_H 
#define __POLYCRYST_H

#include "vector"

#include "Grain.h"
#include "Loading.h"
#include "Params.h"

namespace model
{
	class GrainStructure;

	class Polycrystall
	{
	public:
		
		std::vector<Grain> c;			// Массив элементов поликристалла

		Loading *loading;					// Нагружение
		GrainStructure *structure;			// Зеренная структура поликристалла

		int grainCount;						// Кол-во фрагментов на ребре
		int totalGrainCount;				// Общее кол-во фрагментов

		Tensor D;							// Тензор деформации скорости
		Tensor W;							// Тензор вихря
		Tensor D_in;						// Тензор неупругой части деформации скорости
		Tensor E;							// Тензор деформаций
		Tensor dSgm;						// Тензор скоростей напряжений
		Tensor Sgm;							// Тензор напряжений
		Tensor4 P;							// Усреднённый тензор упругих констант
			
		double stepSavePlot;				// Шаг сохранения графиков
		double stepSavePoleFig;				// Шаг сохранения ПФ
		int stepShowProgress;				// Шаг отображения прогресса
		int stepSaveInternalData;			// Шаг записи отладочных данных
		int periodShowProgress;				// Период обновления процента выполнения

		double strain;						// Интенсивность деформаций
		double stress;						// Интенсивность напряжений

		int fileCount;						// Кол-во отладочных файлов
		std::ofstream *streamDebug;			// Массив файлов для отладки
		std::ofstream *streamInternalVars;	// Массив файлов с кривыми НДС
		std::ofstream *streamDataX;			// Массив файлов с кривыми НДС (X)
		std::ofstream *streamDataY;			// Массив файлов с кривыми НДС (Y)
		std::ofstream *streamDataTest;		// Массив временных (тестовых) файлов

		Polycrystall();
		~Polycrystall();

		void init(int);						// Выделение памяти под зёрна
		void makeGrainStruct();				// Распределение нормалей и фасеток всех фрагментов
		void setParams();					// Распределение параметров фрагментов
		void deformate(Loading*);			// Деформирование поликристалла

		void savePoleFigData();				// Сохранение ПФ
		void saveDebugData();				// Сохранение отладочных данных
		void openAllFiles();				// Открытие всех файлов для записи
		void closeAllFiles();				// Закрытие всех файлов для записи

	private:
		void load(bool unload);				// Нагружение поликристалла
	};

	int get1DPos(int, int, int);			// Возвращает уникальный номер фрагмента в общей структуре
	
	void get3DPos(int, int*, int*, int*);	// Возвращает позицию фрагмента в трёхмерном массиве 

}

#endif __POLYCRYST_H