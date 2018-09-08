#pragma once
#ifndef __ROTATIONS_H 
#define __ROTATIONS_H

#include "Grain.h"

/*
*Ротационные механизмы
*/

namespace model
{
	// Модель стеснённого поворота по Тейлору
	void rotateByTaylor(Grain*);
	// Модель ротаций, связанная с несовместностью сдвигов
	void rotateByTrusov(Grain*);
	// Модель ротационного упрочнения
	void rotationHardening(Grain*);
	
	// Поворот вокруг заданной оси на заданный угол
	void rotate(Grain* f, double dFi, const Vector a);

	/*********************************************************
	********	    Получение полюсных фигур и ССТ	   *******
	*********************************************************/
	//Запись проекций ПФ в файл
	void savePoints(Tensor,	char*, int, int, int);
	//Сохранение полюсной фигуры
	void getPoleFig(Grain*);
	//Запись проекций ССТ в файл
	void saveSSTPoints(Tensor, double, char*, int, int, int);
	//Сохранение ССТ
	void getSST(Grain*);
}

#endif __ROTATIONS_H