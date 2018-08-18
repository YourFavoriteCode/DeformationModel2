﻿#pragma once
#ifndef CORE_VECTOR_H
#define CORE_VECTOR_H

#include "CommonCore.h"

/***********************************************************************
****************                  Вектор               *****************
***********************************************************************/

namespace model
{
	class Vector
	{
	public:
		double c[DIM];					//Компоненты вектора

		double getNorm();				//Возвращает норму вектора
		int normalize();				//Нормализует вектор
		void setZero();					//Обнуляет компоненты вектора
		void set(double,
			double, double);			//Задаёт значения компонент
		double scalMult(Vector v);		//Скалярное произведение векторов
		Vector vectMult(Vector v);		//Векторное произведение векторов

		Vector operator + (Vector);		//Оператор сложения
		Vector operator - (Vector);		//Оператор вычитания
		Vector operator - ();			//Унарный минус
		Vector operator * (double);
		void operator += (Vector);		//Оператор прибавления вектора
		void operator -= (Vector);		//Оператор вычитания вектора
		void operator *= (double);		//Оператор умножения вектора на число
		int operator /= (double);		//Оператор деления вектора на число

		Vector();
		Vector(double, double, double);
		~Vector();

	private:

	};

}

#endif