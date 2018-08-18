#pragma once
#ifndef __LOADING_H 
#define __LOADING_H

#include "Tensor.h"

namespace model
{
	class Loading
	{
	public:

		int cycle;					//Текущий цикл нагружения
		int currentStep;			//Текущий шаг интегрирования

		double tensionComponent;	//Вытягивающая компонента в одноосном нагружении
		double maxStress;			//Значение, до которого разгружать
		double paramLambda;			//Коэффициент в разгрузке
		double additionalStrain;	//Добавочный коэффициент для продолжения циклического нагружения

		Tensor gradV;				//Тензор деформации скорости

		void setLoad(Tensor);

		Loading();
		~Loading();
	private:
	};
}

#endif