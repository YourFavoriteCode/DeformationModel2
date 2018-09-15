#pragma once
#ifndef SLIP_SYSTEM_H
#define SLIP_SYSTEM_H

#include "Tensor.h"

/***********************************************************************
****************         Система скольжения            *****************
***********************************************************************/

namespace model
{

	class SlipSystem
	{
	public:
		Vector n;						// Вектор нормали
		Vector b;						// Вектор Бюргерса
		Tensor o;						// Ориентационный тензор
		double t;						// Действующее касательное напряжение
		double tc;						// Критическое касательное напряжение
		double dgm;						// Скорость сдвига
		double gmm;						// Накопленный сдвиг
		double tbs;						// Обратные касательные напряжения
		
		//Инициализация значений ГЦК/ОЦК
		void Initialize(int, int, int, int, int, int);	
		//Инициализация значений ГПУ
		void Initialize(int, int, int, int,	int, int, int, int, int);

		SlipSystem();
		~SlipSystem();
	};

}

#endif