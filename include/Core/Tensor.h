#pragma once
#ifndef CORE_TENSOR_H
#define CORE_TENSOR_H

#include "CommonCore.h"
#include "Vector.h"

/***********************************************************************
****************                  Тензор               *****************
***********************************************************************/

namespace model
{

	class Tensor
	{
	public:
		double c[DIM][DIM];					//Компоненты тензора

		double getDet();					//Возвращает определитель матрицы компонент
		void set(double, double, double,
			double, double, double,
			double, double, double);		//Задание всех компонент тензора
		void setZero();						//Зануляет компоненты тензора
		void setUnit();						//Делает матрицу компонент тензора единичной
		void transp();						//Транспонирует матрицу компонент тензора
		double doubleScalMult(Tensor);		//Свёртка (двойное скалярное произведение тензоров)
		void rotationMatrix(double, Vector);//Ортогональный тензор поворота вокруг заданной оси на заданный угол
		void getAxisAngle(double*, Vector*);//Возвращает ось и угол поворота, по которым был образован тензор

		Tensor getSymmetryPart();			//Возвращает симметричную часть тензора
		Tensor getAntiSymmetryPart();		//Возвращает антисимметричную часть тензора
		Vector getRow(int);					//Вектор из компонент заданной строки
		Vector getCol(int);					//Вектор из компонент заданного стобца

		double getL(int n);					//Возвращает собственное число с номером n

		Tensor operator + (Tensor);			//Оператор сложения тензоров
		Tensor operator - (Tensor);			//Оператор вычитания тензоров
		Tensor operator - ();				//Унарный минус
		Tensor operator * (Tensor);			//Оператор умножения тензоров

		friend Tensor operator * (Tensor&, double);		//Оператор умножения на число
		friend Tensor operator * (double, Tensor&);		//Коммутативный

		void operator += (Tensor);			//Оператор прибавления тензора
		void operator -= (Tensor);			//Оператор убавления тензора
		void operator *= (Tensor);			//Оператор домножения на тензор
		void operator *= (double);			//Оператор умножения тензора на число
		int operator /= (double);			//Оператор деления тензора на число

		Tensor();
		~Tensor();

	private:

	};
}

#endif