#pragma once
#ifndef __POLYCRYST_H 
#define __POLYCRYST_H

#include "Fragment.h"

namespace model
{
	class Polycrystall
	{
	public:
		Fragment *c;					//Массив элементов поликристалла
		int grainCount;					//Кол-во фрагментов на ребре
		int totalGrainCount;			//Общее кол-во фрагментов

		Tensor D;						//Тензор деформации скорости
		Tensor W;						//Тензор вихря
		Tensor D_in;					//Тензор неупругой части деформации скорости
		Tensor E;						//Тензор деформаций
		Tensor dSgm;					//Тензор скоростей напряжений
		Tensor Sgm;						//Тензор напряжений
		Tensor4 P;						//Усреднённый тензор упругих констант

		int cycle;
		int CURR_STEP;					//Текущий шаг интегрирования
		double PLOT_STEP;				//Шаг сохранения графиков
		double POLUS_STEP;				//Шаг сохранения ПФ
		int PROC_STEP;					//Шаг отображения прогресса
		int DEBUG_STEP;					//Шаг записи отладочных данных
		int proc_period;				//Период обновления процента выполнения

		double strain;					//Интенсивность деформаций
		double stress;					//Интенсивность напряжений

		int fileCount;					//Кол-во отладочных файлов
		std::ofstream *streamDebug;		//Массив файлов для отладки
		std::ofstream *streamInternalVars;	//Массив файлов с кривыми НДС
		std::ofstream *streamDataX;		//Массив файлов с кривыми НДС (X)
		std::ofstream *streamDataY;		//Массив файлов с кривыми НДС (Y)
		std::ofstream *streamDataTest;	//Массив временных (тестовых) файлов

		Polycrystall();
		~Polycrystall();

		double tension_component;	//Вытягивающая компонента в одноосном нагружении
		double final_stress;		//Значение, до которого разгружать
		double lam;					//Коэффициент в разгрузке
		double addition_strain;		//Добавочный коэффициент для продолжения циклического нагружения

		void init(int);					//Выделение памяти под зёрна
		void makeGrainStruct();			//Распределение нормалей и фасеток всех фрагментов
		void setParams();				//Распределение параметров фрагментов
		void deformate();				//Деформирование поликристалла

		void savePoleFigData();			//Сохранение ПФ
		void saveDebugData();			//Сохранение отладочных данных
		void openAllFiles();			//Открытие всех файлов для записи
		void closeAllFiles();			//Закрытие всех файлов для записи

	private:
		void load(bool unload);		//Нагружение поликристалла
	};

	int get1DPos(int, int, int);			//Возвращает уникальный номер фрагмента в общей структуре
	
	void get3DPos(int, int*, int*, int*);	//Возвращает позицию фрагмента в трёхмерном массиве 

}

#endif __POLYCRYST_H