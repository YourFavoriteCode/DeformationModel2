#pragma once
#ifndef __DISTRIB_H 
#define __DISTRIB_H

namespace model
{
	/*
	Различные распределения случайной величины
	Для двухпараметрических: первый параметр - 
	мат.ожидание, второй - дисперсия
	*/

	// Равномерное распределение
	double uniformDistrib(double, double);
	// Нормальное распределение
	double normalDistrib(double m, double d);
	// Логнормальное распределение
	double logNormalDistrib(double m, double d);
	// Показательное распределение
	double expDistrib(double l);
}

#endif __DISTRIB_H