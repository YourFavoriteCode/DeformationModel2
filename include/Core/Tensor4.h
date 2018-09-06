#pragma once
#ifndef CORE_TENSOR4_H
#define CORE_TENSOR4_H

#include "Tensor.h"

/***********************************************************************
**************         Тензор четвёртого ранга            **************
***********************************************************************/

namespace model
{

	class Tensor4
	{
	public:
		double c[DIM][DIM][DIM][DIM];	// Компоненты тензора

		// Обнуление компонент тензора
		void setZero();
		// Симметризация компонент тензора
		void symmetrize();
		// Перевод компонент тензора в ЛСК
		Tensor4 toLsk(Tensor O);
		// Оператор прибавления тензора
		void operator += (Tensor4);
		// Оператор отнимания тензора
		void operator -= (Tensor4);
		// Оператор умножения тензора на число
		void operator *= (double);
		// Оператор деления тензора на число
		int operator /= (double);

		// Двойное скалярное произведение
		Tensor doubleScalMult(Tensor);
		
		Tensor4();
		~Tensor4();

	private:
	};
}

#endif