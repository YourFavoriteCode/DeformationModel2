#pragma once
#ifndef CORE_COMMON_H
#define CORE_COMMON_H

namespace model
{
	/**********************************************************************
	***********				Математические константы			***********
	**********************************************************************/
	const double SQRT3 = 1.732050807568;		//Корень из 3-х
	const double SQRT2 = 1.414213562373;		//Корень из 2-х
	const double SQRT2_3 = 0.816496580927;		//Корень из 2/3
	const double SQRT3_2 = 1.224744871391;		//Корень из 3/2

	const double PI = 3.141592653589;			//Число Пи
	const double PI_2 = 1.570796326794;			//Пи/2
	const double PIx2 = 6.283185307179;			//2Пи

	const double EPS = 1e-10;					//Малая величина

	const int DIM = 3;	//Размерность пространства


	int LeviCivit(int i, int j, int k);	//Псевдо-тензор Леви-Чивита

}

#endif