﻿// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"

#include <cmath>
#include <ctime>
#include <omp.h>
#include <fstream>

#include "Polycrystall.h"
#include "Params.h"
#include "Fragmentation.h"
#include "Distributions.h"
#include "Tension.h"
#include "Rotations.h"
#include "Hardening.h"
#include "Functions.h"

namespace model
{
	Polycrystall::Polycrystall()
	{
		Strain = 0;
		Stress = 0;

		cycle = 0;
		CURR_STEP = 0;
		PROC_STEP = 0;
		POLUS_STEP = 0;
		PLOT_STEP = 0;
		DEBUG_STEP = 0;
		proc_period = 400;
		file_count = 16;

		tension_component = 0;
		final_stress = 1e4;
		lam = 2.5;
		addition_strain = 1e-5;
	}

	Polycrystall::~Polycrystall()
	{
		for (int i = 0; i < fragm_count; i++)
		{
			for (int j = 0; j < fragm_count; j++)
			{
				delete[] C[i][j];
			}
			delete[] C[i];
		}
		delete[] C;
	}

	void Polycrystall::OpenFiles()
	{
		TruncPoleFiles();				//Очистка всех файлов полюсных фигур
		TruncSSTFiles();

		DataXStream = new std::ofstream[20];//Открытие файлов для записи кривых НДС
		DataYStream = new std::ofstream[20];
		Datastream = new std::ofstream[1];

		//Файлы для вывода макро-данных
		if (prms::isSaveMacro && prms::isSaveIntensity)
		{
			DataXStream[0].open("Plot\\macroXint.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[0].open("Plot\\macroYint.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}

		if (prms::isSaveMacro && prms::save11)
		{
			DataXStream[1].open("Plot\\macroX11.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[1].open("Plot\\macroY11.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::isSaveMacro && prms::save12)
		{
			DataXStream[2].open("Plot\\macroX12.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[2].open("Plot\\macroY12.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::isSaveMacro && prms::save13)
		{
			DataXStream[3].open("Plot\\macroX13.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[3].open("Plot\\macroY13.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::isSaveMacro && prms::save21)
		{
			DataXStream[4].open("Plot\\macroX21.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[4].open("Plot\\macroY21.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::isSaveMacro && prms::save22)
		{
			DataXStream[5].open("Plot\\macroX22.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[5].open("Plot\\macroY22.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::isSaveMacro && prms::save23)
		{
			DataXStream[6].open("Plot\\macroX23.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[6].open("Plot\\macroY23.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::isSaveMacro && prms::save31)
		{
			DataXStream[7].open("Plot\\macroX31.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[7].open("Plot\\macroY31.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::isSaveMacro && prms::save32)
		{
			DataXStream[8].open("Plot\\macroX32.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[8].open("Plot\\macroY32.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::isSaveMacro && prms::save33)
		{
			DataXStream[9].open("Plot\\macroX33.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[9].open("Plot\\macroY33.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}

		//Файлы для мезо-данных

		if (prms::isSaveMeso && prms::isSaveIntensity)
		{
			DataXStream[10].open("Plot\\mesoXint.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[10].open("Plot\\mesoYint.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}

		if (prms::isSaveMeso && prms::save11)
		{
			DataXStream[11].open("Plot\\mesoX11.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[11].open("Plot\\mesoY11.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::isSaveMeso && prms::save12)
		{
			DataXStream[12].open("Plot\\mesoX12.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[12].open("Plot\\mesoY12.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::isSaveMeso && prms::save13)
		{
			DataXStream[13].open("Plot\\mesoX13.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[13].open("Plot\\mesoY13.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::isSaveMeso && prms::save21)
		{
			DataXStream[14].open("Plot\\mesoX21.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[14].open("Plot\\mesoY21.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::isSaveMeso && prms::save22)
		{
			DataXStream[15].open("Plot\\mesoX22.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[15].open("Plot\\mesoY22.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::isSaveMeso && prms::save23)
		{
			DataXStream[16].open("Plot\\mesoX23.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[16].open("Plot\\mesoY23.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::isSaveMeso && prms::save31)
		{
			DataXStream[17].open("Plot\\mesoX31.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[17].open("Plot\\mesoY31.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::isSaveMeso && prms::save32)
		{
			DataXStream[18].open("Plot\\mesoX32.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[18].open("Plot\\mesoY32.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::isSaveMeso && prms::save33)
		{
			DataXStream[19].open("Plot\\mesoX33.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[19].open("Plot\\mesoY33.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}

		Datastream[0].open("Plot\\ActiveSS.dat", std::ios_base::out | std::ios_base::trunc | std::ios::binary);

		dbgstream = new std::ofstream[file_count];
		if (prms::debug_period > 0)				//Открытие файлов для отладочных данных
		{
			dbgstream[0].open("DBG\\o.txt", std::ios_base::out | std::ios_base::trunc);
			dbgstream[1].open("DBG\\e.txt", std::ios_base::out | std::ios_base::trunc);
			dbgstream[2].open("DBG\\d.txt", std::ios_base::out | std::ios_base::trunc);
			dbgstream[3].open("DBG\\sgm.txt", std::ios_base::out | std::ios_base::trunc);
			dbgstream[4].open("DBG\\om.txt", std::ios_base::out | std::ios_base::trunc);
			dbgstream[5].open("DBG\\dsgm.txt", std::ios_base::out | std::ios_base::trunc);
			dbgstream[6].open("DBG\\din.txt", std::ios_base::out | std::ios_base::trunc);
			dbgstream[7].open("DBG\\w.txt", std::ios_base::out | std::ios_base::trunc);
			dbgstream[8].open("DBG\\dgamma.txt", std::ios_base::out | std::ios_base::trunc);
			dbgstream[9].open("DBG\\t.txt", std::ios_base::out | std::ios_base::trunc);
			dbgstream[10].open("DBG\\Macro_D.txt", std::ios_base::out | std::ios_base::trunc);
			dbgstream[11].open("DBG\\Macro_Din.txt", std::ios_base::out | std::ios_base::trunc);
			dbgstream[12].open("DBG\\Macro_Sgm.txt", std::ios_base::out | std::ios_base::trunc);
			dbgstream[13].open("DBG\\Macro_dSgm.txt", std::ios_base::out | std::ios_base::trunc);
			dbgstream[14].open("DBG\\Macro_E.txt", std::ios_base::out | std::ios_base::trunc);
			dbgstream[15].open("DBG\\VOL_M.txt", std::ios_base::out | std::ios_base::trunc);
		}

		TestStream = new std::ofstream[6];
		TestStream[0].open("Test0.txt", std::ios_base::out | std::ios_base::trunc);
		TestStream[1].open("Test1.txt", std::ios_base::out | std::ios_base::trunc);
		TestStream[2].open("Test2.txt", std::ios_base::out | std::ios_base::trunc);
		TestStream[3].open("Test3.txt", std::ios_base::out | std::ios_base::trunc);
		TestStream[4].open("Test4.txt", std::ios_base::out | std::ios_base::trunc);
		TestStream[5].open("Test5.txt", std::ios_base::out | std::ios_base::trunc);
	}

	void Polycrystall::CloseFiles()
	{
		for (int i = 0; i < 20; i++)
		{
			DataXStream[i].close();
			DataYStream[i].close();
		}
		Datastream[0].close();
		for (int i = 0; i < 6; i++)
		{
			TestStream[i].close();
		}

		if (prms::debug_period > 0)
		{
			for (int i = 0; i < file_count; i++)
			{
				dbgstream[i].close();
			}
		}
	}

	void Polycrystall::Init(int count)
	{
		fragm_count = count;
		total_fragm_count = (int)pow(count, 3);

		C = new Fragment**[count];		//Выделение памяти под массив
		for (int i = 0; i < count; i++)
		{
			C[i] = new Fragment*[count];
			for (int j = 0; j < count; j++)
			{
				C[i][j] = new Fragment[count];
			}
		}
	}

	void Polycrystall::setParams()
	{
		for (int q1 = 0; q1 < fragm_count; q1++)
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{
					//Задание материала 
					int another_material;//Примесная фаза
					another_material = (prms::materialType == 1) ? 0 : 1;

					int a = (int)(((double)rand() / RAND_MAX) * 100);//На всё воля божья
					if (a <= prms::material_purity)
					{
						C[q1][q2][q3].setMaterialParams(prms::materialType);
					}
					else
					{
						C[q1][q2][q3].setMaterialParams(another_material);
					}

					C[q1][q2][q3].rot_Mc = prms::rotationParamMc;	//Раздача начальных критических моментов
					C[q1][q2][q3].rot_A = prms::rotationParamA;	//и параметров модели ротаций
					C[q1][q2][q3].rot_H = prms::rotationParamH;
					C[q1][q2][q3].rot_L = prms::rotationParamL;
					C[q1][q2][q3].position = get1DPos(q1, q2, q3);//Получение порядкового номера фрагмента

					if (prms::orientationType == 0)
					{
						//Углами Эйлера
						if (prms::randomOrientations)//Получение ориентационного тензора (случайный равномерный закон распределения)
						{
							double a = ((double)rand() / RAND_MAX) * (PI);
							double g = ((double)rand() / RAND_MAX) * (PI);
							double y1 = ((double)rand() / RAND_MAX);
							double y2 = ((double)rand() / RAND_MAX);
							double cb = y1 > 0.5 ? y2 : -y2;
							C[q1][q2][q3].Orientate(a, g, cb);
						}
						else//Получение ориентационного тензора (КСК=ЛСК)
						{
							C[q1][q2][q3].o.setUnit();
						}
					}
					else if (prms::orientationType == 1)
					{
						//С помощью оси и угла
						if (prms::randomOrientations)//Получение ориентационного тензора (случайный равномерный закон распределения)
						{
							double fi = ((double)rand() / RAND_MAX) * (PI);
							double psi = ((double)rand() / RAND_MAX) * (PIx2);
							double y1 = ((double)rand() / RAND_MAX);
							double y2 = ((double)rand() / RAND_MAX);
							//Равномерное распределение косинуса угла
							double cf = y1 > 0.5 ? y2 : -y2;
							//Переход из сферической в ортогональную СК
							double x = sin(fi)*cos(psi);
							double y = sin(fi)*sin(psi);
							double z = cos(fi);
							//Запись оси
							Vector axis;
							axis.set(x, y, z);
							axis.Normalize();
							C[q1][q2][q3].OrientateAxis(cf, axis);
						}
						else//Получение ориентационного тензора (КСК=ЛСК)
						{
							C[q1][q2][q3].o.setUnit();
						}
					}
					else if (prms::orientationType == 2)
					{
						//С помощью кватерниона
						if (prms::randomOrientations)//Получение ориентационного тензора (случайный равномерный закон распределения)
						{
							double fi = ((double)rand() / RAND_MAX) * (PI);
							double psi = ((double)rand() / RAND_MAX) * (PIx2);
							//Равномерное распределение косинуса
							double y1 = ((double)rand() / RAND_MAX);
							double y2 = ((double)rand() / RAND_MAX);
							double w = y1 > 0.5 ? y2 : -y2;
							double buf = sqrt(1.0 - w * w);
							double x = sin(fi)*cos(psi)*buf;
							double y = sin(fi)*sin(psi)*buf;
							double z = cos(fi)*buf;
							C[q1][q2][q3].OrientateQuater(w, x, y, z);
						}
						else//Получение ориентационного тензора (КСК=ЛСК)
						{
							C[q1][q2][q3].o.setUnit();
						}
					}

					//Задание размеров фрагментов
					switch (prms::fragm_size_law)
					{
					case 0://Равномерное
					{
						C[q1][q2][q3].size = UniformDistrib(prms::fragm_size_m, prms::fragm_size_dsp);
						break;
					}
					case 1://Нормальное
					{
						C[q1][q2][q3].size = NormalDistrib(prms::fragm_size_m, prms::fragm_size_dsp);
						break;
					}
					case 2://Логнормальное
					{
						C[q1][q2][q3].size = LogNormalDistrib(prms::fragm_size_m, prms::fragm_size_dsp);
						break;
					}
					case 3://Показательное
					{
						C[q1][q2][q3].size = ExpDistrib(prms::fragm_size_m);//Только один параметр
						break;
					}
					}
					C[q1][q2][q3].volume = pow(C[q1][q2][q3].size, 3);	//Объём фрагмента

					//Выделение памяти под массивы, необходимые для работы с окружением
					C[q1][q2][q3].surrounds = new Fragment[prms::surroundCount];
					C[q1][q2][q3].normals = new Vector[prms::surroundCount];
					C[q1][q2][q3].contact = new int[prms::surroundCount];

					for (int h = 0; h < prms::surroundCount; h++)
					{
						C[q1][q2][q3].contact[h] = -1;		//Изначально контакт не задан
					}
				}
			}
		}
	}

	void Polycrystall::MakeStruct()
	{
		for (int q1 = 0; q1 < fragm_count; q1++)
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{

					for (int h = 0; h < prms::surroundCount; h++)
					{
						//Если контакт уже был задан - пропускаем
						if (C[q1][q2][q3].contact[h] != -1) continue;
						//Определяем, граничат ли фрагменты
						//Первые 6, т.е. боковые грани, граничат всегда
						double a = h < 6 ? 1 : ((double)rand() / RAND_MAX);//На всё воля божья
						if (a < 0.5)
						{
							//Контакта нет - тоже пропускаем
							C[q1][q2][q3].contact[h] = 0;
							continue;
						}

						int qq1 = q1, qq2 = q2, qq3 = q3, y;
						//qq1, qq2, qq3 - координаты зерна соседа
						//y - номер нормали в соседнем зерне в направлении данного зерна
						double fi = ((double)rand() / RAND_MAX) * (PI / 12);//Случайный угол отклонения нормали
						//TODO: предвычислить наиболее распространенные слагаемые для удобства чтения
						switch (h)
						{
						case 0://Вверх
						{
							C[q1][q2][q3].normals[h].set(-sin(fi), sin(fi) / cos(fi), 1 / cos(fi));
							qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
							y = 5;
							break;
						}
						case 1://От нас
						{
							C[q1][q2][q3].normals[h].set(-1 / cos(fi), sin(fi), sin(fi) / cos(fi));
							qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
							y = 3;
							break;
						}
						case 2://Вправо
						{
							C[q1][q2][q3].normals[h].set(sin(fi) / cos(fi), 1 / cos(fi), -sin(fi));
							qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
							y = 4;
							break;
						}
						case 3://На нас
						{
							C[q1][q2][q3].normals[h].set(1 / cos(fi), -sin(fi), sin(fi) / cos(fi));
							qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
							y = 1;
							break;
						}
						case 4://Влево
						{
							C[q1][q2][q3].normals[h].set(sin(fi), -1 / cos(fi), -sin(fi) / cos(fi));
							qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
							y = 2;
							break;
						}
						case 5://Вниз
						{
							C[q1][q2][q3].normals[h].set(sin(fi) / cos(fi), sin(fi), -1 / cos(fi));
							qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
							y = 0;
							break;
						}
						//Далее идут уже необязательные соседи
						/**************           Рёбра куба          ***************************/
						case 6://Лево от нас
						{
							C[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * cos(fi)*cos(PI_2), -cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
							qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
							qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
							y = 9;
							break;
						}
						case 7://Лево на нас
						{
							C[q1][q2][q3].normals[h].set(cos(fi + PI_2) * cos(fi)*cos(PI_2), -cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
							qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
							qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
							y = 8;
							break;
						}
						case 8://право от нас
						{
							C[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
							qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
							qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
							y = 7;
							break;
						}
						case 9://Право на нас
						{
							C[q1][q2][q3].normals[h].set(cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
							qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
							qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
							y = 6;
							break;
						}
						case 10://Верх лево
						{
							C[q1][q2][q3].normals[h].set(cos(fi + PI_2) * sin(fi), -cos(fi + PI_2) * cos(fi), sin(fi + PI_2)*cos(fi));
							qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
							qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
							y = 15;
							break;
						}
						case 11://Верх право
						{
							C[q1][q2][q3].normals[h].set(cos(fi + PI_2) * sin(fi), cos(fi + PI_2) * cos(fi), sin(fi + PI_2)*cos(fi));
							qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
							qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
							y = 14;
							break;
						}
						case 12://Верх на нас
						{
							C[q1][q2][q3].normals[h].set(cos(fi + PI_2) * cos(fi), cos(fi + PI_2) * sin(fi), sin(fi + PI_2)*cos(fi));
							qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
							qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
							y = 17;
							break;
						}
						case 13://Верх от нас
						{
							C[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * cos(fi), cos(fi + PI_2) * sin(fi), sin(fi + PI_2)*cos(fi));
							qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
							qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
							y = 16;
							break;
						}
						case 14://Низ лево
						{
							C[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * sin(fi), -cos(fi + PI_2) * cos(fi), -sin(fi + PI_2)*cos(fi));
							qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
							qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
							y = 11;
							break;
						}
						case 15://Низ право
						{
							C[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * sin(fi), cos(fi + PI_2) * cos(fi), -sin(fi + PI_2)*cos(fi));
							qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
							qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
							y = 10;
							break;
						}
						case 16://Низ на нас
						{
							C[q1][q2][q3].normals[h].set(cos(fi + PI_2) * cos(fi), -cos(fi + PI_2) * sin(fi), -sin(fi + PI_2)*cos(fi));
							qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
							qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
							y = 13;
							break;
						}
						case 17://Низ от нас
						{
							C[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * cos(fi), -cos(fi + PI_2) * sin(fi), -sin(fi + PI_2)*cos(fi));
							qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
							qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
							y = 12;
							break;
						}
						/**************      Вершины     *****************/
						case 18://верх лево от нас
						{
							C[q1][q2][q3].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
							qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
							qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
							qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
							y = 25;
							break;
						}
						case 19://верх лево на нас
						{
							C[q1][q2][q3].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
							qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
							qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
							qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
							y = 24;
							break;
						}
						case 20://верх право от нас
						{
							C[q1][q2][q3].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
							qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
							qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
							qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
							y = 23;
							break;
						}
						case 21://верх право на нас
						{
							C[q1][q2][q3].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
							qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
							qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
							qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
							y = 22;
							break;
						}
						case 22://низ лево от нас
						{
							C[q1][q2][q3].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
							qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
							qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
							qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
							y = 21;
							break;
						}
						case 23://низ лево на нас
						{
							C[q1][q2][q3].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
							qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
							qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
							qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
							y = 20;
							break;
						}
						case 24://низ право от нас
						{
							C[q1][q2][q3].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));

							qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
							qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
							qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
							y = 19;
							break;
						}
						case 25://низ право на нас
						{
							C[q1][q2][q3].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
							qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
							qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
							qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
							y = 18;
							break;
						}
						}

						C[q1][q2][q3].surrounds[h] = C[qq1][qq2][qq3];//Здравствуй, сосед!
						C[qq1][qq2][qq3].surrounds[y] = C[q1][q2][q3];//Приятно познакомиться!
						C[q1][q2][q3].normals[h].Normalize();

						for (int i = 0; i < DIM; i++)
						{
							C[qq1][qq2][qq3].normals[y].C[i] = -C[q1][q2][q3].normals[h].C[i];//Поделись нормалью
						}

						if (h < 6) C[q1][q2][q3].contact[h] = 1;		//Контакт на грани октаэдра
						else if (h < 14) C[q1][q2][q3].contact[h] = 3;	//Контакт на вершине октаэдра
						else C[q1][q2][q3].contact[h] = 2;				//Контакт на ребре
					}
					if (prms::surroundCount > 6)	//Уменьшение объёма из-за отсечений
					{
						double a = C[q1][q2][q3].size * 0.1;			//Длина срезанной части вдоль ребра
						double vol_edge = a*a*C[q1][q2][q3].size / 2.0;	//Объём, срезанный рёбрами
						double vol_vertex = a*a*a / SQRT3;				//Объём, срезанный вершинами
						int cut_edge = 0;		//Кол-во срезанных рёбер
						int cut_vertex = 0;		//Кол-во срезанных вершин
						for (int h = 6; h < prms::surroundCount; h++)
						{
							if (C[q1][q2][q3].contact[h] != 0)
							{
								if (h < 14) cut_vertex++;
								else cut_edge++;
							}
						}
						C[q1][q2][q3].volume -= (cut_edge*vol_edge + cut_vertex*vol_vertex);//Вычитание
					}

				}
			}
		}
	}

	void Polycrystall::SavePoleFig()
	{
		for (int q1 = 0; q1 < fragm_count; q1++)
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{
					GetPoleFig(&C[q1][q2][q3]);
					if (prms::SST_SAVING) GetSST(&C[q1][q2][q3]);
				}
			}
		}
	}

	void Polycrystall::SaveDbgInfo()
	{
		for (int i = 0; i < file_count; i++)//Визуальное разделение шагов 
		{
			dbgstream[i] << "#########################      STEP " << CURR_STEP << "      #########################" << std::endl << std::endl;
		}

		//Запись тензоров каждого из зерен или фрагментов
		for (int q1 = 0; q1 < fragm_count; q1++)
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{
					WriteDebugInfo(dbgstream[0], C[q1][q2][q3].o.C);
					WriteDebugInfo(dbgstream[1], C[q1][q2][q3].e.C);
					WriteDebugInfo(dbgstream[2], C[q1][q2][q3].d.C);
					WriteDebugInfo(dbgstream[3], C[q1][q2][q3].sgm.C);
					WriteDebugInfo(dbgstream[4], C[q1][q2][q3].om.C);
					WriteDebugInfo(dbgstream[5], C[q1][q2][q3].dsgm.C);
					WriteDebugInfo(dbgstream[6], C[q1][q2][q3].d_in.C);
					WriteDebugInfo(dbgstream[7], C[q1][q2][q3].w.C);
					for (int f = 0; f < C[q1][q2][q3].SS_count; f++)
					{
						dbgstream[8] << C[q1][q2][q3].SS[f].dgm << " ";
					}
					dbgstream[8] << std::endl << std::endl;
					for (int f = 0; f < C[q1][q2][q3].SS_count; f++)
					{
						dbgstream[9] << C[q1][q2][q3].SS[f].t << " ";
					}
					dbgstream[9] << std::endl << std::endl;
					dbgstream[15] << C[q1][q2][q3].moment.C[0] << " " << C[q1][q2][q3].moment.C[1] << " " << C[q1][q2][q3].moment.C[2] << std::endl;
				}
			}
		}
		//Запись тензоров представительного объема
		WriteDebugInfo(dbgstream[10], D.C);
		WriteDebugInfo(dbgstream[11], D_in.C);
		WriteDebugInfo(dbgstream[12], Sgm.C);
		WriteDebugInfo(dbgstream[13], dSgm.C);
		WriteDebugInfo(dbgstream[14], E.C);
	}

	void Polycrystall::Load(bool unload)
	{
		E += D*prms::dt;
		/*Параметр unload включает разгрузку представительного объёма*/
		if (prms::trueUniaxial || unload)	//Одноосное растяжение
		{
			//Осреднение
			P.setZero();
			D_in.setZero();
			for (int q1 = 0; q1 < fragm_count; q1++)
			{
				for (int q2 = 0; q2 < fragm_count; q2++)
				{
					for (int q3 = 0; q3 < fragm_count; q3++)
					{
						D_in += C[q1][q2][q3].d_in;
						P += C[q1][q2][q3].p.ToLSK(C[q1][q2][q3].o);
					}
				}
			}

			D_in /= (total_fragm_count);
			P /= (total_fragm_count);

			//Симметризация тензора упругих констант
			P.Symmetrize();

			D = !unload ? TensionStrainCalc(P, D_in, D.C[0][0]) : UnloadingStrainCalc(P, D_in, Sgm, lam);

			Strain = SQRT2_3*sqrt(E.doubleScalMult(E));//Вычисление интенсивности деформаций

			dSgm = TensionStressCalc(P, D_in, D);
			//dSgm *= prms::dt;				//Приращение напряжений на шаге
			Sgm += dSgm*prms::dt;
			Stress = SQRT3_2*sqrt(Sgm.doubleScalMult(Sgm));//Вычисление интенсивности напряжений
		}
		else
		{
			Stress = 0;		//Вычисление интенсивностей осреднением
			Strain = 0;
			for (int q1 = 0; q1 < fragm_count; q1++)
			{
				for (int q2 = 0; q2 < fragm_count; q2++)
				{
					for (int q3 = 0; q3 < fragm_count; q3++)
					{
						Strain += C[q1][q2][q3].strain;
						Stress += C[q1][q2][q3].stress;
					}
				}
			}
			Strain /= total_fragm_count;
			Stress /= total_fragm_count;
		}

#pragma omp parallel for
		//Часть, которую можно паралелить
		//Здесь необходимо гарантировать защиту данных каждого фрагмента
		//от перезаписи другими фрагментами
		for (int q1 = 0; q1 < fragm_count; q1++)
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{

					/**************************************************
					************       Переходим в КСК       **********
					**************************************************/

					Tensor O = C[q1][q2][q3].o;
					Tensor OT = O;
					OT.Transp();
					C[q1][q2][q3].d = O*D*OT;//Гипотеза Фойгта
					C[q1][q2][q3].w = O*W*OT /*- C[q1][q2][q3].om*/;//Расширенная

					C[q1][q2][q3].sgm = O*C[q1][q2][q3].sgm*OT;
					C[q1][q2][q3].d_in = O*C[q1][q2][q3].d_in*OT;


					/***************************************************
					***********       Пересчитываем НДС      ***********
					***************************************************/

					C[q1][q2][q3].NDScalc();

					if (prms::usingHardeningBase)			//Базовое упрочнение
					{
						Base_hardening(&C[q1][q2][q3]);
					}

					if (prms::ROTATIONS_TAYLOR)		//Ротации по Тейлору
					{
						Taylor_rotations(&C[q1][q2][q3]);
					}

					if (prms::ROTATIONS_TRUSOV && prms::ROTATIONS_HARDENING)	//Ротационное упрочнение
					{
						Rotation_hardening(&C[q1][q2][q3]);
					}

					if (prms::usingHardeningBound)	//Зернограничное упрочнение
					{
						Boundary_hardening(&C[q1][q2][q3]);
					}

					if (prms::ROTATIONS_TRUSOV)		//Ротации по Трусову
					{
						Trusov_rotations(&C[q1][q2][q3]);
					}
					/**************************************************
					************       Переходим в ЛСК       **********
					**************************************************/

					C[q1][q2][q3].sgm = OT*C[q1][q2][q3].sgm*O;
					C[q1][q2][q3].d_in = OT*C[q1][q2][q3].d_in*O;
					C[q1][q2][q3].iter++;
				}
			}
		}


		if (!prms::trueUniaxial && !unload)		//Этот блок нужен исключительно для работы с энергией!
		{
			Sgm.setZero();
			D_in.setZero();
			for (int q1 = 0; q1 < fragm_count; q1++)
			{
				for (int q2 = 0; q2 < fragm_count; q2++)
				{
					for (int q3 = 0; q3 < fragm_count; q3++)
					{
						Sgm += C[q1][q2][q3].sgm;
						D_in += C[q1][q2][q3].d_in;
					}
				}
			}
			Sgm /= total_fragm_count;
			D_in /= total_fragm_count;
		}


		/************************************************************
		***********	        Прогресс выполнения 	      ***********
		************************************************************/

		double progress;
		if (!unload)
		{
			progress = prms::trueUniaxial ? fabs(E.C[0][0]) : Strain;
			progress = progress / prms::maxStrainIntencity * 100.0;

			if (!(prms::cycle_count == 1 || cycle == 0))	//Для многоцикловых нагружений
			{
				progress /= 2.0;
				/************************************************
				*Этот код позволяет корректно отображать прогресс
				*выполнения во время циклических нагружений
				*************************************************/
				if (E.C[0][0] > 0)
				{
					if (Sgm.C[0][0] > 0)
					{
						progress += 50.0;
					}
					else
					{
						progress = 50.0 - progress;
					}
				}
				else
				{
					if (Sgm.C[0][0] > 0)
					{
						progress = 50.0 - progress;
					}
					else
					{
						progress += 50.0;
					}
				}
			}
		}
		else
		{
			progress = final_stress / fabs(Sgm.C[0][0]) * 100.0;	//Индикация прогресса при разгрузке
		}

		int period = unload ? proc_period / 40 : proc_period;

		if (PROC_STEP == period)
		{
			PROC_STEP = 0;
			//Курсор двигается на 6 символов влево, записывается новое значение с точностью 5.2 и знак %
			printf("\b\b\b\b\b\b%0*.*f%%", 5, 2, progress);
		}

		/************************************************************
		***********	    Запись данных для графиков НДС    ***********
		************************************************************/

		if ((progress - PLOT_STEP > prms::plot_period || unload) && prms::plot_period > 0)
		{
			BoundsAnalize();		//Подсчет доли большеугловых границ

			if (prms::isSaveMacro)	//Запись компонент тензоров макроуровня
			{
				if (prms::isSaveIntensity)
				{
					DataXStream[0].write((char *)&Strain, sizeof(double));
					DataYStream[0].write((char *)&Stress, sizeof(double));
				}
				if (prms::save11)
				{
					DataXStream[1].write((char *)&E.C[0][0], sizeof(double));
					DataYStream[1].write((char *)&Sgm.C[0][0], sizeof(double));
				}
				if (prms::save12)
				{
					DataXStream[2].write((char *)&E.C[0][1], sizeof(double));
					DataYStream[2].write((char *)&Sgm.C[0][1], sizeof(double));
				}
				if (prms::save13)
				{
					DataXStream[3].write((char *)&E.C[0][2], sizeof(double));
					DataYStream[3].write((char *)&Sgm.C[0][2], sizeof(double));
				}
				if (prms::save21)
				{
					DataXStream[4].write((char *)&E.C[1][0], sizeof(double));
					DataYStream[4].write((char *)&Sgm.C[1][0], sizeof(double));
				}
				if (prms::save22)
				{
					DataXStream[5].write((char *)&E.C[1][1], sizeof(double));
					DataYStream[5].write((char *)&Sgm.C[1][1], sizeof(double));
				}
				if (prms::save23)
				{
					DataXStream[6].write((char *)&E.C[1][2], sizeof(double));
					DataYStream[6].write((char *)&Sgm.C[1][2], sizeof(double));
				}
				if (prms::save31)
				{
					DataXStream[7].write((char *)&E.C[2][0], sizeof(double));
					DataYStream[7].write((char *)&Sgm.C[2][0], sizeof(double));
				}
				if (prms::save32)
				{
					DataXStream[8].write((char *)&E.C[2][1], sizeof(double));
					DataYStream[8].write((char *)&Sgm.C[2][1], sizeof(double));
				}
				if (prms::save33)
				{
					DataXStream[9].write((char *)&E.C[2][2], sizeof(double));
					DataYStream[9].write((char *)&Sgm.C[2][2], sizeof(double));
				}
			}

			if (prms::isSaveMeso)	//Запись компонент тензоров мезоуровня
			{
				for (int q1 = 0; q1 < fragm_count; q1++)
				{
					for (int q2 = 0; q2 < fragm_count; q2++)
					{
						for (int q3 = 0; q3 < fragm_count; q3++)
						{

							if (prms::isSaveIntensity)
							{
								DataXStream[10].write((char *)&C[q1][q2][q3].strain, sizeof(double));
								DataYStream[10].write((char *)&C[q1][q2][q3].stress, sizeof(double));
							}
							if (prms::save11)
							{
								DataXStream[11].write((char *)&C[q1][q2][q3].e.C[0][0], sizeof(double));
								DataYStream[11].write((char *)&C[q1][q2][q3].sgm.C[0][0], sizeof(double));
							}
							if (prms::save12)
							{
								DataXStream[12].write((char *)&C[q1][q2][q3].e.C[0][1], sizeof(double));
								DataYStream[12].write((char *)&C[q1][q2][q3].sgm.C[0][1], sizeof(double));
							}
							if (prms::save13)
							{
								DataXStream[13].write((char *)&C[q1][q2][q3].e.C[0][2], sizeof(double));
								DataYStream[13].write((char *)&C[q1][q2][q3].sgm.C[0][2], sizeof(double));
							}
							if (prms::save21)
							{
								DataXStream[14].write((char *)&C[q1][q2][q3].e.C[1][0], sizeof(double));
								DataYStream[14].write((char *)&C[q1][q2][q3].sgm.C[1][0], sizeof(double));
							}
							if (prms::save22)
							{
								DataXStream[15].write((char *)&C[q1][q2][q3].e.C[1][1], sizeof(double));
								DataYStream[15].write((char *)&C[q1][q2][q3].sgm.C[1][1], sizeof(double));
							}
							if (prms::save23)
							{
								DataXStream[16].write((char *)&C[q1][q2][q3].e.C[1][2], sizeof(double));
								DataYStream[16].write((char *)&C[q1][q2][q3].sgm.C[1][2], sizeof(double));
							}
							if (prms::save31)
							{
								DataXStream[17].write((char *)&C[q1][q2][q3].e.C[2][0], sizeof(double));
								DataYStream[17].write((char *)&C[q1][q2][q3].sgm.C[2][0], sizeof(double));
							}
							if (prms::save32)
							{
								DataXStream[18].write((char *)&C[q1][q2][q3].e.C[2][1], sizeof(double));
								DataYStream[18].write((char *)&C[q1][q2][q3].sgm.C[2][1], sizeof(double));
							}
							if (prms::save33)
							{
								DataXStream[19].write((char *)&C[q1][q2][q3].e.C[2][2], sizeof(double));
								DataYStream[19].write((char *)&C[q1][q2][q3].sgm.C[2][2], sizeof(double));
							}

						}
					}
				}
			}

			double ActiveSysCount = 0;			//Среднее кол-во активных систем скольжения на шаге
			double RotEnergy = 0;				//Энергия ротаций на шаге
			double RotSpeed = 0;				//Средняя скорость вращения на шаге
			int RotCount = 0;					//Кол-во вращающихся фрагментов
			double norma = 0;
			double Mc = 0;
			double dmc = 0;
			double angle = 0;
			double H = 0;
			for (int q1 = 0; q1 < fragm_count; q1++)
			{
				for (int q2 = 0; q2 < fragm_count; q2++)
				{
					for (int q3 = 0; q3 < fragm_count; q3++)
					{
						for (int i = 0; i < C[q1][q2][q3].SS_count; i++)
						{
							if (C[q1][q2][q3].SS[i].dgm > EPS) ActiveSysCount++;//Подсчёт активных СС
						}

						if (C[q1][q2][q3].isRotate) RotCount++;		//Подсчёт вращающихся решёток
						RotEnergy += C[q1][q2][q3].rot_energy;		//Суммирование энергий вращения
						RotSpeed += C[q1][q2][q3].rot_speed;		//Суммирование скоростей вращения
						norma += C[q1][q2][q3].norm;
						Mc += C[q1][q2][q3].rot_Mc;
						dmc += C[q1][q2][q3].dmc;
						angle += C[q1][q2][q3].sum_angle;
						H += C[q1][q2][q3].rot_H;
					}
				}
			}
			H /= total_fragm_count;
			angle /= total_fragm_count;
			norma /= total_fragm_count;
			Mc /= total_fragm_count;
			dmc /= total_fragm_count;
			ActiveSysCount /= total_fragm_count;
			if (prms::isSaveActiveSS) Datastream[0].write((char *)&ActiveSysCount, sizeof ActiveSysCount);//Запись кол-ва активных СС
			if (RotCount != 0)
			{
				RotSpeed /= RotCount;
			}
			else RotSpeed = 0;

			/*******************************************************
			********* 		   Работа с энергией          **********
			*******************************************************/

			//Полная энергия деформирования - сумма элементарных энергий на каждом шаге
			//Элементарная энергия - свёртка напряжений с приращением деформации
			//Энергия ротаций - момент*приращение угла

			Tensor dE = D;
			dE *= prms::dt;			//Приращение деформации на шаге

			double StepEnergy = Sgm.doubleScalMult(D);	//Полная энергия на шаге
			double StepEnergy_in = Sgm.doubleScalMult(D_in);	//Полная энергия на шаге
			TestStream[0] << RotCount << std::endl;
			TestStream[1] << RotSpeed << std::endl;
			TestStream[2] << RotEnergy << std::endl;
			TestStream[3] << StepEnergy << std::endl;
			TestStream[4] << StepEnergy_in << std::endl;
			TestStream[5] << norma << std::endl;

			PLOT_STEP = progress;
		}

		/************************************************************
		***********	      Сохранение полюсных фигур	      ***********
		************************************************************/
		if (progress - POLUS_STEP > prms::polus_period && prms::polus_period > 0)
		{
			SavePoleFig();
			POLUS_STEP = progress;
		}

		/************************************************************
		***********	       Запись пошаговых данных	      ***********
		************************************************************/
		if (CURR_STEP >= prms::DEBUG_START && CURR_STEP <= prms::DEBUG_STOP && DEBUG_STEP == prms::debug_period)
		{
			DEBUG_STEP = 0;
			SaveDbgInfo();
		}
		CURR_STEP++;
		PROC_STEP++;
		if (CURR_STEP >= prms::DEBUG_START && CURR_STEP <= prms::DEBUG_STOP) DEBUG_STEP++;

	}

	void Polycrystall::Deformate()
		/****************************************
		*Функция деформирования представительного
		*объема поликристала.
		*****************************************/
	{
		omp_set_num_threads(prms::thread_count);	//Кол-во используемых потоков

		if (prms::trueUniaxial)
		{
			/***************************************************
			*Выбор растягивающей компоненты тензора D.
			*Относительно неё на каждом шаге будет решаться СЛАУ
			***************************************************/
			tension_component = D.C[0][0];
		}
		for (cycle = 0; cycle < prms::cycle_count; cycle++)
		{
			unsigned long t1, t2;

			PLOT_STEP = 0;		//Обнуление счетчиков для периодического вывода данных в файлы
			POLUS_STEP = 0;
			PROC_STEP = 0;
			if (prms::cycle_count > 1)
			{
				printf("\n Loading #%d", cycle + 1);
				t1 = clock();		//Начальная отсечка времени
			}
			printf("\n        00.00%");
			double* counter;
			counter = !prms::trueUniaxial ? &Strain : &E.C[0][0];
			//Для циклических нагружений цикл ведется по значению растягивающей компоненты
			//а для остальных - по значению интенсивности тензора деформации

			while (fabs(*counter) < prms::maxStrainIntencity)	//Цикл по деформациям
			{
				Load(false);
				//if (prms::FRAGMENTATION) GrainRotate();
			}
			if (prms::cycle_count > 1)
			{
				t2 = clock();		//Конечная отсечка времени
				printf("\b\b\b\b\b\bDone in %g sec", (t2 - t1) / 1000.0);
			}
			//		printf("\n START 2! \n");
			//	D.set(0, 0.003, 0, 0.003, 0, 0, 0, 0, 0);//Простой сдвиг
			//	W.setZero();
			//	D.set(-0.003, 0.003, 0, -0.003, 0.003, 0, 0, 0, 0);//РКУП
			/*	D.set(0.003, 0, 0, 0, -0.0015, 0, 0, 0, -0.0015);//Одноосье
				prms::strain_max += prms::strain_max;

				while (fabs(*counter) < prms::strain_max)
				{
				Load(false);
				}
				*/

			if (prms::withUnloading)	//Упругая разгрузка
			{
				printf("\n Unloading #%d", cycle + 1);
				printf("\n        00.00%");
				t1 = clock();		//Начальная отсечка времени
				while (fabs(Sgm.C[0][0]) > final_stress) //Цикл по напряжениям
				{
					Load(true);
				}
				t2 = clock();		//Конечная отсечка времени
				printf("\b\b\b\b\b\bDone in %g sec", (t2 - t1) / 1000.0);
			}

			if (prms::cycle_count > 1)	//Цикоическое знакопеременное нагружение
			{
				D.C[0][0] = pow(-1, cycle + 1) * tension_component;	//Меняем знак растягивающей компоненты
				prms::maxStrainIntencity += prms::maxStrainIntencity * addition_strain;	//Повышаем предел интенсивности
			}

			if (prms::FRAGMENTATION) Fragmentate();

		}
	}
}