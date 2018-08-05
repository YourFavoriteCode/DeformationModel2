#pragma once
#ifndef CORE_FUNCTIONS_H
#define CORE_FUNCTIONS_H

#include "Vector.h"
#include "Tensor.h"

namespace model
{
	int LeviCivit(int i, int j, int k);	//Псевдо-тензор Леви-Чивита

	Tensor VectMult(Vector, Tensor);	//Векторное произведение вектора на тензор
	Tensor VectMult(Tensor, Vector);	//Векторное произведение тензора на вектор
	Vector ScalMult(Vector, Tensor);	//Скалярное произведение вектора на тензор
	Vector ScalMult(Tensor, Vector);	//Скалярное произведение тензора на вектор

	Tensor Transp(Tensor);				//Возвращает транспонированный аргумент
}

#endif