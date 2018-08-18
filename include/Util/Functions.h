#pragma once
#ifndef __FUNCTIONS_H 
#define __FUNCTIONS_H

#include "windows.h"
#include <fstream>

namespace model
{
	/*****************************************************
	*********	 Функции для работы с файлами	 *********
	******************************************************/
	bool isDirectoryExists(LPCWSTR);						// Проверка на существование директории
	void writeDebugInfo(std::ofstream&, double [3][3]);		// Запись данных в файл
	void truncPoleFigFiles();								// Очистка файлов полюсных фигур
	void truncSSTFiles();									// Очистка файлов стандартных треугольников

	bool isNormalDouble(double);							// Проверка на валидность числа
}

#endif __FUNCTIONS_H