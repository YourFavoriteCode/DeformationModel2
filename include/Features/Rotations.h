#pragma once
#ifndef __ROTATIONS_H 
#define __ROTATIONS_H

#include "Grain.h"

/*
*Ротационные механизмы
*/

namespace model
{
	
	void Taylor_rotations(Grain*);		//Модель стеснённого поворота по Тейлору
	void Trusov_rotations(Grain*);		//Модель, связанная с несовместностью сдвигов
	void Rotation_hardening(Grain*);		//Модель ротационного упрочнения
	void Rotate(Grain* f, double dFi, const Vector a);	//Поворот вокруг заданной оси на заданный угол
	/*********************************************************
	********	    Получение полюсных фигур и ССТ	   *******
	*********************************************************/

	void SavePoints(Tensor,	char*,
		int, int, int);			//Запись проекций ПФ в файл
	void GetPoleFig(Grain*);				//Сохранение полюсной фигуры
	
	void SaveSSTPoints(Tensor&, float,		//Запись проекций ССТ в файл
		char*, int, int, int);
	void GetSST(Grain*);					//Сохранение ССТ
}

#endif __ROTATIONS_H