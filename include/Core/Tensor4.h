#pragma once
#ifndef CORE_TENSOR4_H
#define CORE_TENSOR4_H

#include "CommonCore.h"
#include "Tensor.h"

/***********************************************************************
**************         Тензор четвёртого ранга            **************
***********************************************************************/

namespace model
{

	class Tensor4
	{
	public:
		double C[DIM][DIM][DIM][DIM];	//Компоненты тензора

		void setZero();					//Обнуление компонент тензора
		void Symmetrize();				//Симметризация компонент тензора

		Tensor4 ToLSK(Tensor O);		//Перевод компонент тензора в ЛСК

		void operator += (Tensor4);		//Оператор прибавления тензора
		void operator -= (Tensor4);		//Оператор отнимания тензора
		void operator *= (double);		//Оператор умножения тензора на число
		int operator /= (double);		//Оператор деления тензора на число

		Tensor4();
		~Tensor4();
	private:
	};
}

#endif