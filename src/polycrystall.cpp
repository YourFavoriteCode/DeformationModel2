// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"

#include <cmath>
#include <ctime>
#include <omp.h>
#include <fstream>

#include "Polycrystall.h"
#include "Params.h"
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
		delete[] C;
	}

	void Polycrystall::OpenFiles()
	{
		truncPoleFigFiles();				//Очистка всех файлов полюсных фигур
		truncSSTFiles();

		DataXStream = new std::ofstream[20];//Открытие файлов для записи кривых НДС
		DataYStream = new std::ofstream[20];
		Datastream = new std::ofstream[1];

		//Файлы для вывода макро-данных
		if (prms::saveMacro && prms::saveIntensity)
		{
			DataXStream[0].open("Plot\\macroXint.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[0].open("Plot\\macroYint.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}

		if (prms::saveMacro && prms::save11)
		{
			DataXStream[1].open("Plot\\macroX11.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[1].open("Plot\\macroY11.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMacro && prms::save12)
		{
			DataXStream[2].open("Plot\\macroX12.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[2].open("Plot\\macroY12.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMacro && prms::save13)
		{
			DataXStream[3].open("Plot\\macroX13.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[3].open("Plot\\macroY13.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMacro && prms::save21)
		{
			DataXStream[4].open("Plot\\macroX21.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[4].open("Plot\\macroY21.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMacro && prms::save22)
		{
			DataXStream[5].open("Plot\\macroX22.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[5].open("Plot\\macroY22.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMacro && prms::save23)
		{
			DataXStream[6].open("Plot\\macroX23.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[6].open("Plot\\macroY23.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMacro && prms::save31)
		{
			DataXStream[7].open("Plot\\macroX31.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[7].open("Plot\\macroY31.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMacro && prms::save32)
		{
			DataXStream[8].open("Plot\\macroX32.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[8].open("Plot\\macroY32.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMacro && prms::save33)
		{
			DataXStream[9].open("Plot\\macroX33.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[9].open("Plot\\macroY33.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}

		//Файлы для мезо-данных

		if (prms::saveMeso && prms::saveIntensity)
		{
			DataXStream[10].open("Plot\\mesoXint.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[10].open("Plot\\mesoYint.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}

		if (prms::saveMeso && prms::save11)
		{
			DataXStream[11].open("Plot\\mesoX11.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[11].open("Plot\\mesoY11.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMeso && prms::save12)
		{
			DataXStream[12].open("Plot\\mesoX12.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[12].open("Plot\\mesoY12.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMeso && prms::save13)
		{
			DataXStream[13].open("Plot\\mesoX13.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[13].open("Plot\\mesoY13.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMeso && prms::save21)
		{
			DataXStream[14].open("Plot\\mesoX21.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[14].open("Plot\\mesoY21.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMeso && prms::save22)
		{
			DataXStream[15].open("Plot\\mesoX22.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[15].open("Plot\\mesoY22.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMeso && prms::save23)
		{
			DataXStream[16].open("Plot\\mesoX23.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[16].open("Plot\\mesoY23.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMeso && prms::save31)
		{
			DataXStream[17].open("Plot\\mesoX31.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[17].open("Plot\\mesoY31.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMeso && prms::save32)
		{
			DataXStream[18].open("Plot\\mesoX32.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[18].open("Plot\\mesoY32.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMeso && prms::save33)
		{
			DataXStream[19].open("Plot\\mesoX33.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			DataYStream[19].open("Plot\\mesoY33.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}

		Datastream[0].open("Plot\\ActiveSS.dat", std::ios_base::out | std::ios_base::trunc | std::ios::binary);

		dbgstream = new std::ofstream[file_count];
		if (prms::saveVariablesPeriodStep > 0)				//Открытие файлов для отладочных данных
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

		if (prms::saveVariablesPeriodStep > 0)
		{
			for (int i = 0; i < file_count; i++)
			{
				dbgstream[i].close();
			}
		}
	}

	void Polycrystall::Init(int count)
	{
		grainCount = count;
		totalGrainCount = (int)pow(count, 3);

		C = new Fragment[totalGrainCount];		//Выделение памяти под массив
	}

	void Polycrystall::setParams()
	{
		for (int q = 0; q < totalGrainCount; q++)
		{

			//Задание материала 
			int another_material;//Примесная фаза
			another_material = (prms::materialType == 1) ? 0 : 1;

			int a = (int)(((double)rand() / RAND_MAX) * 100);//На всё воля божья
			if (a <= prms::mainPhasePercent)
			{
				C[q].setMaterialParams(prms::materialType);
			}
			else
			{
				C[q].setMaterialParams(another_material);
			}

			C[q].rot_Mc = prms::rotationParamMc;	//Раздача начальных критических моментов
			C[q].rot_A = prms::rotationParamA;	//и параметров модели ротаций
			C[q].rot_H = prms::rotationParamH;
			C[q].rot_L = prms::rotationParamL;
			C[q].position = q;//Получение порядкового номера фрагмента

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
					C[q].Orientate(a, g, cb);
				}
				else//Получение ориентационного тензора (КСК=ЛСК)
				{
					C[q].o.setUnit();
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
					axis.normalize();
					C[q].OrientateAxis(cf, axis);
				}
				else//Получение ориентационного тензора (КСК=ЛСК)
				{
					C[q].o.setUnit();
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
					C[q].OrientateQuater(w, x, y, z);
				}
				else//Получение ориентационного тензора (КСК=ЛСК)
				{
					C[q].o.setUnit();
				}
			}

			//Задание размеров фрагментов
			switch (prms::grainSizeDistribLaw)
			{
			case prms::DISTRIB_UNIFORM:
			{
				C[q].size = UniformDistrib(prms::grainSizeDistribM, prms::grainSizeDistribD);
				break;
			}
			case prms::DISTRIB_NORMAL:
			{
				C[q].size = NormalDistrib(prms::grainSizeDistribM, prms::grainSizeDistribD);
				break;
			}
			case prms::DISTRIB_LOGNORMAL:
			{
				C[q].size = LogNormalDistrib(prms::grainSizeDistribM, prms::grainSizeDistribD);
				break;
			}
			case prms::DISTRIB_EXPONENT:
			{
				C[q].size = ExpDistrib(prms::grainSizeDistribM);//Только один параметр
				break;
			}
			}
			C[q].volume = pow(C[q].size, 3);	//Объём фрагмента

			//Выделение памяти под массивы, необходимые для работы с окружением
			C[q].surrounds = new Fragment[prms::grainSurroundCount];
			C[q].normals = new Vector[prms::grainSurroundCount];
			C[q].contact = new int[prms::grainSurroundCount];

			for (int h = 0; h < prms::grainSurroundCount; h++)
			{
				C[q].contact[h] = -1;		//Изначально контакт не задан
			}

		}
	}

	void Polycrystall::MakeStruct()
	{
		for (int q = 0; q < totalGrainCount; q++)
		{
			int q1, q2, q3, y;
			get3DPos(q, &q1, &q2, &q3);
			for (int h = 0; h < prms::grainSurroundCount; h++)
			{
				//Если контакт уже был задан - пропускаем
				if (C[q].contact[h] != -1) continue;
				//Определяем, граничат ли фрагменты
				//Первые 6, т.е. боковые грани, граничат всегда
				double a = h < 6 ? 1 : ((double)rand() / RAND_MAX);//На всё воля божья
				if (a < 0.5)
				{
					//Контакта нет - тоже пропускаем
					C[q].contact[h] = 0;
					continue;
				}
				int qq1 = q1, qq2 = q2, qq3 = q3;
				//qq1, qq2, qq3 - координаты зерна соседа
				//y - номер нормали в соседнем зерне в направлении данного зерна
				double fi = ((double)rand() / RAND_MAX) * (PI / 12);//Случайный угол отклонения нормали
				//TODO: предвычислить наиболее распространенные слагаемые для удобства чтения
				switch (h)
				{
				case 0://Вверх
				{
					C[q].normals[h].set(-sin(fi), sin(fi) / cos(fi), 1 / cos(fi));
					qq3 = q3 == grainCount - 1 ? 0 : q3 + 1;
					y = 5;
					break;
				}
				case 1://От нас
				{
					C[q].normals[h].set(-1 / cos(fi), sin(fi), sin(fi) / cos(fi));
					qq1 = q1 == 0 ? grainCount - 1 : q1 - 1;
					y = 3;
					break;
				}
				case 2://Вправо
				{
					C[q].normals[h].set(sin(fi) / cos(fi), 1 / cos(fi), -sin(fi));
					qq2 = q2 == grainCount - 1 ? 0 : q2 + 1;
					y = 4;
					break;
				}
				case 3://На нас
				{
					C[q].normals[h].set(1 / cos(fi), -sin(fi), sin(fi) / cos(fi));
					qq1 = q1 == grainCount - 1 ? 0 : q1 + 1;
					y = 1;
					break;
				}
				case 4://Влево
				{
					C[q].normals[h].set(sin(fi), -1 / cos(fi), -sin(fi) / cos(fi));
					qq2 = q2 == 0 ? grainCount - 1 : q2 - 1;
					y = 2;
					break;
				}
				case 5://Вниз
				{
					C[q].normals[h].set(sin(fi) / cos(fi), sin(fi), -1 / cos(fi));
					qq3 = q3 == 0 ? grainCount - 1 : q3 - 1;
					y = 0;
					break;
				}
				//Далее идут уже необязательные соседи
				/**************           Рёбра куба          ***************************/
				case 6://Лево от нас
				{
					C[q].normals[h].set(-cos(fi + PI_2) * cos(fi)*cos(PI_2), -cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
					qq1 = q1 == 0 ? grainCount - 1 : q1 - 1;
					qq2 = q2 == 0 ? grainCount - 1 : q2 - 1;
					y = 9;
					break;
				}
				case 7://Лево на нас
				{
					C[q].normals[h].set(cos(fi + PI_2) * cos(fi)*cos(PI_2), -cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
					qq1 = q1 == grainCount - 1 ? 0 : q1 + 1;
					qq2 = q2 == 0 ? grainCount - 1 : q2 - 1;
					y = 8;
					break;
				}
				case 8://право от нас
				{
					C[q].normals[h].set(-cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
					qq1 = q1 == 0 ? grainCount - 1 : q1 - 1;
					qq2 = q2 == grainCount - 1 ? 0 : q2 + 1;
					y = 7;
					break;
				}
				case 9://Право на нас
				{
					C[q].normals[h].set(cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
					qq2 = q2 == grainCount - 1 ? 0 : q2 + 1;
					qq1 = q1 == grainCount - 1 ? 0 : q1 + 1;
					y = 6;
					break;
				}
				case 10://Верх лево
				{
					C[q].normals[h].set(cos(fi + PI_2) * sin(fi), -cos(fi + PI_2) * cos(fi), sin(fi + PI_2)*cos(fi));
					qq3 = q3 == grainCount - 1 ? 0 : q3 + 1;
					qq2 = q2 == 0 ? grainCount - 1 : q2 - 1;
					y = 15;
					break;
				}
				case 11://Верх право
				{
					C[q].normals[h].set(cos(fi + PI_2) * sin(fi), cos(fi + PI_2) * cos(fi), sin(fi + PI_2)*cos(fi));
					qq3 = q3 == grainCount - 1 ? 0 : q3 + 1;
					qq2 = q2 == grainCount - 1 ? 0 : q2 + 1;
					y = 14;
					break;
				}
				case 12://Верх на нас
				{
					C[q].normals[h].set(cos(fi + PI_2) * cos(fi), cos(fi + PI_2) * sin(fi), sin(fi + PI_2)*cos(fi));
					qq3 = q3 == grainCount - 1 ? 0 : q3 + 1;
					qq1 = q1 == grainCount - 1 ? 0 : q1 + 1;
					y = 17;
					break;
				}
				case 13://Верх от нас
				{
					C[q].normals[h].set(-cos(fi + PI_2) * cos(fi), cos(fi + PI_2) * sin(fi), sin(fi + PI_2)*cos(fi));
					qq3 = q3 == grainCount - 1 ? 0 : q3 + 1;
					qq1 = q1 == 0 ? grainCount - 1 : q1 - 1;
					y = 16;
					break;
				}
				case 14://Низ лево
				{
					C[q].normals[h].set(-cos(fi + PI_2) * sin(fi), -cos(fi + PI_2) * cos(fi), -sin(fi + PI_2)*cos(fi));
					qq2 = q2 == 0 ? grainCount - 1 : q2 - 1;
					qq3 = q3 == 0 ? grainCount - 1 : q3 - 1;
					y = 11;
					break;
				}
				case 15://Низ право
				{
					C[q].normals[h].set(-cos(fi + PI_2) * sin(fi), cos(fi + PI_2) * cos(fi), -sin(fi + PI_2)*cos(fi));
					qq2 = q2 == grainCount - 1 ? 0 : q2 + 1;
					qq3 = q3 == 0 ? grainCount - 1 : q3 - 1;
					y = 10;
					break;
				}
				case 16://Низ на нас
				{
					C[q].normals[h].set(cos(fi + PI_2) * cos(fi), -cos(fi + PI_2) * sin(fi), -sin(fi + PI_2)*cos(fi));
					qq1 = q1 == grainCount - 1 ? 0 : q1 + 1;
					qq3 = q3 == 0 ? grainCount - 1 : q3 - 1;
					y = 13;
					break;
				}
				case 17://Низ от нас
				{
					C[q].normals[h].set(-cos(fi + PI_2) * cos(fi), -cos(fi + PI_2) * sin(fi), -sin(fi + PI_2)*cos(fi));
					qq1 = q1 == 0 ? grainCount - 1 : q1 - 1;
					qq3 = q3 == 0 ? grainCount - 1 : q3 - 1;
					y = 12;
					break;
				}
				/**************      Вершины     *****************/
				case 18://верх лево от нас
				{
					C[q].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
					qq1 = q1 == 0 ? grainCount - 1 : q1 - 1;
					qq2 = q2 == 0 ? grainCount - 1 : q2 - 1;
					qq3 = q3 == grainCount - 1 ? 0 : q3 + 1;
					y = 25;
					break;
				}
				case 19://верх лево на нас
				{
					C[q].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
					qq1 = q1 == grainCount - 1 ? 0 : q1 + 1;
					qq2 = q2 == 0 ? grainCount - 1 : q2 - 1;
					qq3 = q3 == grainCount - 1 ? 0 : q3 + 1;
					y = 24;
					break;
				}
				case 20://верх право от нас
				{
					C[q].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
					qq1 = q1 == 0 ? grainCount - 1 : q1 - 1;
					qq2 = q2 == grainCount - 1 ? 0 : q2 + 1;
					qq3 = q3 == grainCount - 1 ? 0 : q3 + 1;
					y = 23;
					break;
				}
				case 21://верх право на нас
				{
					C[q].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
					qq1 = q1 == grainCount - 1 ? 0 : q1 + 1;
					qq2 = q2 == grainCount - 1 ? 0 : q2 + 1;
					qq3 = q3 == grainCount - 1 ? 0 : q3 + 1;
					y = 22;
					break;
				}
				case 22://низ лево от нас
				{
					C[q].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
					qq1 = q1 == 0 ? grainCount - 1 : q1 - 1;
					qq2 = q2 == 0 ? grainCount - 1 : q2 - 1;
					qq3 = q3 == 0 ? grainCount - 1 : q3 - 1;
					y = 21;
					break;
				}
				case 23://низ лево на нас
				{
					C[q].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
					qq1 = q1 == grainCount - 1 ? 0 : q1 + 1;
					qq2 = q2 == 0 ? grainCount - 1 : q2 - 1;
					qq3 = q3 == 0 ? grainCount - 1 : q3 - 1;
					y = 20;
					break;
				}
				case 24://низ право от нас
				{
					C[q].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));

					qq1 = q1 == 0 ? grainCount - 1 : q1 - 1;
					qq2 = q2 == grainCount - 1 ? 0 : q2 + 1;
					qq3 = q3 == 0 ? grainCount - 1 : q3 - 1;
					y = 19;
					break;
				}
				case 25://низ право на нас
				{
					C[q].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
					qq1 = q1 == grainCount - 1 ? 0 : q1 + 1;
					qq2 = q2 == grainCount - 1 ? 0 : q2 + 1;
					qq3 = q3 == 0 ? grainCount - 1 : q3 - 1;
					y = 18;
					break;
				}
				}
				int index = get1DPos(qq1, qq2, qq3);
				C[q].surrounds[h] = C[index];//Здравствуй, сосед!
				C[index].surrounds[y] = C[q];//Приятно познакомиться!
				C[q].normals[h].normalize();

				for (int i = 0; i < DIM; i++)
				{
					C[index].normals[y].c[i] = -C[q].normals[h].c[i];//Поделись нормалью
				}

				if (h < 6) C[q].contact[h] = 1;		//Контакт на грани октаэдра
				else if (h < 14) C[q].contact[h] = 3;	//Контакт на вершине октаэдра
				else C[q].contact[h] = 2;				//Контакт на ребре
			}
			if (prms::grainSurroundCount > 6)	//Уменьшение объёма из-за отсечений
			{
				double a = C[q].size * 0.1;			//Длина срезанной части вдоль ребра
				double vol_edge = a * a*C[q].size / 2.0;	//Объём, срезанный рёбрами
				double vol_vertex = a * a*a / SQRT3;				//Объём, срезанный вершинами
				int cut_edge = 0;		//Кол-во срезанных рёбер
				int cut_vertex = 0;		//Кол-во срезанных вершин
				for (int h = 6; h < prms::grainSurroundCount; h++)
				{
					if (C[q].contact[h] != 0)
					{
						if (h < 14) cut_vertex++;
						else cut_edge++;
					}
				}
				C[q].volume -= (cut_edge*vol_edge + cut_vertex * vol_vertex);//Вычитание
			}


		}
	}

	void Polycrystall::SavePoleFig()
	{
		for (int q = 0; q < totalGrainCount; q++)
		{
			GetPoleFig(&C[q]);
			if (prms::usingStandardTriangleSaving) GetSST(&C[q]);

		}
	}

	void Polycrystall::SaveDbgInfo()
	{
		for (int i = 0; i < file_count; i++)//Визуальное разделение шагов 
		{
			dbgstream[i] << "#########################      STEP " << CURR_STEP << "      #########################" << std::endl << std::endl;
		}

		//Запись тензоров каждого из зерен или фрагментов
		for (int q = 0; q < totalGrainCount; q++)
		{

			writeDebugInfo(dbgstream[0], C[q].o.c);
			writeDebugInfo(dbgstream[1], C[q].e.c);
			writeDebugInfo(dbgstream[2], C[q].d.c);
			writeDebugInfo(dbgstream[3], C[q].sgm.c);
			writeDebugInfo(dbgstream[4], C[q].om.c);
			writeDebugInfo(dbgstream[5], C[q].dsgm.c);
			writeDebugInfo(dbgstream[6], C[q].d_in.c);
			writeDebugInfo(dbgstream[7], C[q].w.c);
			for (int f = 0; f < C[q].SS_count; f++)
			{
				dbgstream[8] << C[q].SS[f].dgm << " ";
			}
			dbgstream[8] << std::endl << std::endl;
			for (int f = 0; f < C[q].SS_count; f++)
			{
				dbgstream[9] << C[q].SS[f].t << " ";
			}
			dbgstream[9] << std::endl << std::endl;
			dbgstream[15] << C[q].moment.c[0] << " " << C[q].moment.c[1] << " " << C[q].moment.c[2] << std::endl;

		}
		//Запись тензоров представительного объема
		writeDebugInfo(dbgstream[10], D.c);
		writeDebugInfo(dbgstream[11], D_in.c);
		writeDebugInfo(dbgstream[12], Sgm.c);
		writeDebugInfo(dbgstream[13], dSgm.c);
		writeDebugInfo(dbgstream[14], E.c);
	}

	void Polycrystall::Load(bool unload)
	{
		E += D * prms::dt;
		/*Параметр unload включает разгрузку представительного объёма*/
		if (prms::trueUniaxial || unload)	//Одноосное растяжение
		{
			//Осреднение
			P.setZero();
			D_in.setZero();
			for (int q = 0; q < totalGrainCount; q++)
			{
				D_in += C[q].d_in;
				P += C[q].p.ToLSK(C[q].o);
			}

			D_in /= (totalGrainCount);
			P /= (totalGrainCount);

			//Симметризация тензора упругих констант
			P.Symmetrize();

			D = !unload ? TensionStrainCalc(P, D_in, D.c[0][0]) : UnloadingStrainCalc(P, D_in, Sgm, lam);

			Strain = SQRT2_3 * sqrt(E.doubleScalMult(E));//Вычисление интенсивности деформаций

			dSgm = TensionStressCalc(P, D_in, D);
			//dSgm *= prms::dt;				//Приращение напряжений на шаге
			Sgm += dSgm * prms::dt;
			Stress = SQRT3_2 * sqrt(Sgm.doubleScalMult(Sgm));//Вычисление интенсивности напряжений
		}
		else
		{
			Stress = 0;		//Вычисление интенсивностей осреднением
			Strain = 0;
			for (int q = 0; q < totalGrainCount; q++)
			{
				Strain += C[q].strain;
				Stress += C[q].stress;

			}
			Strain /= totalGrainCount;
			Stress /= totalGrainCount;
		}

#pragma omp parallel for
		//Часть, которую можно паралелить
		//Здесь необходимо гарантировать защиту данных каждого фрагмента
		//от перезаписи другими фрагментами
		for (int q = 0; q < totalGrainCount; q++)
		{

			/**************************************************
			************       Переходим в КСК       **********
			**************************************************/

			Tensor O = C[q].o;
			Tensor OT = O;
			OT.transp();
			C[q].d = O * D*OT;//Гипотеза Фойгта
			C[q].w = O * W*OT /*- C[q].om*/;//Расширенная

			C[q].sgm = O * C[q].sgm*OT;
			C[q].d_in = O * C[q].d_in*OT;


			/***************************************************
			***********       Пересчитываем НДС      ***********
			***************************************************/

			C[q].NDScalc();

			if (prms::usingHardeningBase)			//Базовое упрочнение
			{
				Base_hardening(&C[q]);
			}

			if (prms::usingRotationsTaylor)		//Ротации по Тейлору
			{
				Taylor_rotations(&C[q]);
			}

			if (prms::usingRotationsTrusov && prms::usingRotationsHardening)	//Ротационное упрочнение
			{
				Rotation_hardening(&C[q]);
			}

			if (prms::usingHardeningBound)	//Зернограничное упрочнение
			{
				Boundary_hardening(&C[q]);
			}

			if (prms::usingRotationsTrusov)		//Ротации по Трусову
			{
				Trusov_rotations(&C[q]);
			}
			/**************************************************
			************       Переходим в ЛСК       **********
			**************************************************/

			C[q].sgm = OT * C[q].sgm*O;
			C[q].d_in = OT * C[q].d_in*O;
			C[q].iter++;

		}


		if (!prms::trueUniaxial && !unload)		//Этот блок нужен исключительно для работы с энергией!
		{
			Sgm.setZero();
			D_in.setZero();
			for (int q = 0; q < totalGrainCount; q++)
			{

				Sgm += C[q].sgm;
				D_in += C[q].d_in;

			}
			Sgm /= totalGrainCount;
			D_in /= totalGrainCount;
		}


		/************************************************************
		***********	        Прогресс выполнения 	      ***********
		************************************************************/

		double progress;
		if (!unload)
		{
			progress = prms::trueUniaxial ? fabs(E.c[0][0]) : Strain;
			progress = progress / prms::maxStrainIntencity * 100.0;

			if (!(prms::loadCycleCount == 1 || cycle == 0))	//Для многоцикловых нагружений
			{
				progress /= 2.0;
				/************************************************
				*Этот код позволяет корректно отображать прогресс
				*выполнения во время циклических нагружений
				*************************************************/
				if (E.c[0][0] > 0)
				{
					if (Sgm.c[0][0] > 0)
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
					if (Sgm.c[0][0] > 0)
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
			progress = final_stress / fabs(Sgm.c[0][0]) * 100.0;	//Индикация прогресса при разгрузке
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

		if ((progress - PLOT_STEP > prms::periodSavePlot || unload) && prms::periodSavePlot > 0)
		{

			if (prms::saveMacro)	//Запись компонент тензоров макроуровня
			{
				if (prms::saveIntensity)
				{
					DataXStream[0].write((char *)&Strain, sizeof(double));
					DataYStream[0].write((char *)&Stress, sizeof(double));
				}
				if (prms::save11)
				{
					DataXStream[1].write((char *)&E.c[0][0], sizeof(double));
					DataYStream[1].write((char *)&Sgm.c[0][0], sizeof(double));
				}
				if (prms::save12)
				{
					DataXStream[2].write((char *)&E.c[0][1], sizeof(double));
					DataYStream[2].write((char *)&Sgm.c[0][1], sizeof(double));
				}
				if (prms::save13)
				{
					DataXStream[3].write((char *)&E.c[0][2], sizeof(double));
					DataYStream[3].write((char *)&Sgm.c[0][2], sizeof(double));
				}
				if (prms::save21)
				{
					DataXStream[4].write((char *)&E.c[1][0], sizeof(double));
					DataYStream[4].write((char *)&Sgm.c[1][0], sizeof(double));
				}
				if (prms::save22)
				{
					DataXStream[5].write((char *)&E.c[1][1], sizeof(double));
					DataYStream[5].write((char *)&Sgm.c[1][1], sizeof(double));
				}
				if (prms::save23)
				{
					DataXStream[6].write((char *)&E.c[1][2], sizeof(double));
					DataYStream[6].write((char *)&Sgm.c[1][2], sizeof(double));
				}
				if (prms::save31)
				{
					DataXStream[7].write((char *)&E.c[2][0], sizeof(double));
					DataYStream[7].write((char *)&Sgm.c[2][0], sizeof(double));
				}
				if (prms::save32)
				{
					DataXStream[8].write((char *)&E.c[2][1], sizeof(double));
					DataYStream[8].write((char *)&Sgm.c[2][1], sizeof(double));
				}
				if (prms::save33)
				{
					DataXStream[9].write((char *)&E.c[2][2], sizeof(double));
					DataYStream[9].write((char *)&Sgm.c[2][2], sizeof(double));
				}
			}

			if (prms::saveMeso)	//Запись компонент тензоров мезоуровня
			{
				for (int q = 0; q < totalGrainCount; q++)
				{

					if (prms::saveIntensity)
					{
						DataXStream[10].write((char *)&C[q].strain, sizeof(double));
						DataYStream[10].write((char *)&C[q].stress, sizeof(double));
					}
					if (prms::save11)
					{
						DataXStream[11].write((char *)&C[q].e.c[0][0], sizeof(double));
						DataYStream[11].write((char *)&C[q].sgm.c[0][0], sizeof(double));
					}
					if (prms::save12)
					{
						DataXStream[12].write((char *)&C[q].e.c[0][1], sizeof(double));
						DataYStream[12].write((char *)&C[q].sgm.c[0][1], sizeof(double));
					}
					if (prms::save13)
					{
						DataXStream[13].write((char *)&C[q].e.c[0][2], sizeof(double));
						DataYStream[13].write((char *)&C[q].sgm.c[0][2], sizeof(double));
					}
					if (prms::save21)
					{
						DataXStream[14].write((char *)&C[q].e.c[1][0], sizeof(double));
						DataYStream[14].write((char *)&C[q].sgm.c[1][0], sizeof(double));
					}
					if (prms::save22)
					{
						DataXStream[15].write((char *)&C[q].e.c[1][1], sizeof(double));
						DataYStream[15].write((char *)&C[q].sgm.c[1][1], sizeof(double));
					}
					if (prms::save23)
					{
						DataXStream[16].write((char *)&C[q].e.c[1][2], sizeof(double));
						DataYStream[16].write((char *)&C[q].sgm.c[1][2], sizeof(double));
					}
					if (prms::save31)
					{
						DataXStream[17].write((char *)&C[q].e.c[2][0], sizeof(double));
						DataYStream[17].write((char *)&C[q].sgm.c[2][0], sizeof(double));
					}
					if (prms::save32)
					{
						DataXStream[18].write((char *)&C[q].e.c[2][1], sizeof(double));
						DataYStream[18].write((char *)&C[q].sgm.c[2][1], sizeof(double));
					}
					if (prms::save33)
					{
						DataXStream[19].write((char *)&C[q].e.c[2][2], sizeof(double));
						DataYStream[19].write((char *)&C[q].sgm.c[2][2], sizeof(double));
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
			for (int q = 0; q < totalGrainCount; q++)
			{
				for (int i = 0; i < C[q].SS_count; i++)
				{
					if (C[q].SS[i].dgm > EPS) ActiveSysCount++;//Подсчёт активных СС
				}

				if (C[q].isRotate) RotCount++;		//Подсчёт вращающихся решёток
				RotEnergy += C[q].rot_energy;		//Суммирование энергий вращения
				RotSpeed += C[q].rot_speed;		//Суммирование скоростей вращения
				norma += C[q].norm;
				Mc += C[q].rot_Mc;
				dmc += C[q].dmc;
				angle += C[q].sum_angle;
				H += C[q].rot_H;
			}
			H /= totalGrainCount;
			angle /= totalGrainCount;
			norma /= totalGrainCount;
			Mc /= totalGrainCount;
			dmc /= totalGrainCount;
			ActiveSysCount /= totalGrainCount;
			if (prms::saveActiveSS) Datastream[0].write((char *)&ActiveSysCount, sizeof ActiveSysCount);//Запись кол-ва активных СС
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
		if (progress - POLUS_STEP > prms::periodSavePolus && prms::periodSavePolus > 0)
		{
			SavePoleFig();
			POLUS_STEP = progress;
		}

		/************************************************************
		***********	       Запись пошаговых данных	      ***********
		************************************************************/
		if (CURR_STEP >= prms::saveVariablesStartStep && CURR_STEP <= prms::saveVariablesStopStep && DEBUG_STEP == prms::saveVariablesPeriodStep)
		{
			DEBUG_STEP = 0;
			SaveDbgInfo();
		}
		CURR_STEP++;
		PROC_STEP++;
		if (CURR_STEP >= prms::saveVariablesStartStep && CURR_STEP <= prms::saveVariablesStopStep) DEBUG_STEP++;

	}

	void Polycrystall::Deformate()
		/****************************************
		*Функция деформирования представительного
		*объема поликристала.
		*****************************************/
	{
		omp_set_num_threads(prms::ompThreadCount);	//Кол-во используемых потоков

		if (prms::trueUniaxial)
		{
			/***************************************************
			*Выбор растягивающей компоненты тензора D.
			*Относительно неё на каждом шаге будет решаться СЛАУ
			***************************************************/
			tension_component = D.c[0][0];
		}
		for (cycle = 0; cycle < prms::loadCycleCount; cycle++)
		{
			unsigned long t1, t2;

			PLOT_STEP = 0;		//Обнуление счетчиков для периодического вывода данных в файлы
			POLUS_STEP = 0;
			PROC_STEP = 0;
			if (prms::loadCycleCount > 1)
			{
				printf("\n Loading #%d", cycle + 1);
				t1 = clock();		//Начальная отсечка времени
			}
			printf("\n        00.00%");
			double* counter;
			counter = !prms::trueUniaxial ? &Strain : &E.c[0][0];
			//Для циклических нагружений цикл ведется по значению растягивающей компоненты
			//а для остальных - по значению интенсивности тензора деформации

			while (fabs(*counter) < prms::maxStrainIntencity)	//Цикл по деформациям
			{
				Load(false);
				//if (prms::FRAGMENTATION) GrainRotate();
			}
			if (prms::loadCycleCount > 1)
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
				while (fabs(Sgm.c[0][0]) > final_stress) //Цикл по напряжениям
				{
					Load(true);
				}
				t2 = clock();		//Конечная отсечка времени
				printf("\b\b\b\b\b\bDone in %g sec", (t2 - t1) / 1000.0);
			}

			if (prms::loadCycleCount > 1)	//Цикоическое знакопеременное нагружение
			{
				D.c[0][0] = pow(-1, cycle + 1) * tension_component;	//Меняем знак растягивающей компоненты
				prms::maxStrainIntencity += prms::maxStrainIntencity * addition_strain;	//Повышаем предел интенсивности
			}

		}
	}

	int get1DPos(int q1, int q2, int q3)
	{
		//По трём координатам в объёме поликристалла возвращает уникальный номер фрагмента
		int res = q1 * prms::grainCountLinear*prms::grainCountLinear + q2 * prms::grainCountLinear + q3;
		return res;
	}

	void get3DPos(int pos, int* q1, int* q2, int* q3)
	{
		//Восстанавливает пространственные координаты поликристалла по уникальному номеру
		int C2d = prms::grainCountLinear*prms::grainCountLinear;
		int C3d = C2d * prms::grainCountLinear;

		int qq1 = pos / C2d;
		int qq2 = (pos - qq1 * C2d) / prms::grainCountLinear;
		int qq3 = (pos - qq1 * C2d) % prms::grainCountLinear;

		*q1 = qq1;
		*q2 = qq2;
		*q3 = qq3;
	}

}