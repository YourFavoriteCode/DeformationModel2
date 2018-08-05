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
	bool isDirectoryExists(LPCWSTR);						//Проверка на существование директории
	void WriteDebugInfo(std::ofstream&, double [3][3]);		//Запись данных в файл
	void TruncPoleFiles();									//Очистка файлов полюсных фигур
	void TruncSSTFiles();									//Очистка файлов стандартных треугольников

}

#endif __FUNCTIONS_H