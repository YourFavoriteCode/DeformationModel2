// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"

#include <fstream>
#include <ctime>
#include <Windows.h>

#include "Params.h"
#include "Functions.h"
#include "Polycrystall.h"

using namespace model;

int _tmain(int argc, _TCHAR* argv[])
{
	if (argc == 1) return 1;	//Программа закроется, если вызвана без аргументов
	/****************************************************************
	*********	  Создание несуществующих директорий		*********
	****************************************************************/
	if (!isDirectoryExists(L"Plot"))
	{
		CreateDirectory((LPCTSTR)"Plot", NULL);//Для графиков
	}
	if (!isDirectoryExists(L"Polus"))
	{
		CreateDirectory((LPCTSTR)"Polus", NULL);//Для полюсных фигур
	}
	if (!isDirectoryExists(L"DBG"))
	{
		CreateDirectory((LPCTSTR)"DBG", NULL);//Для отладочных данных
	}

	/*****************************************************************
	********     Интерфейс ввода/вывода параметров модели     ********
	*****************************************************************/
	printf(" Build date %s, %s\n", __DATE__, __TIME__);
	char* param_file = new char[256];
	//wcstombs(param_file, argv[1], 256);//Получили имя файла с параметрами
	param_file = argv[1];
	printf(" Parameters file: %s\n", param_file);
	if (prms::ReadParams(param_file) == 1) printf(" Error in file!\n");		//Считали параметры из файла
	delete[] param_file;//Больше не нужен
	printf(" ==========================================\n");
	const int total_fragm_count = (int)pow(prms::fragm_count, 3);	//Общее кол-во фрагментов
	printf(" Fragments count: %d\n", total_fragm_count);
	printf(" Max. strain: %g\n", prms::maxStrainIntencity);
	if (prms::cycle_count != 1)
	{
		printf(" Number of cycles: %d\n", prms::cycle_count);
	}
	printf(" Integration step: %g\n", prms::dt);
	if (prms::ROTATIONS_HARDENING)
	{
		printf(" Using rotation hardening model\n");
		printf("       K1: %g\n", prms::rotationParamHardK1);
		printf("       K2: %g\n", prms::rotationParamHardK2);
	}
	if (prms::ROTATIONS_TAYLOR)
	{
		printf(" Using Taylor's rotation model\n");
	}
	if (prms::ROTATIONS_TRUSOV)
	{
		printf(" Using Trusov's rotation model\n");
		printf("       A: %g\n", prms::rotationParamA);
		printf("       H: %g\n", prms::rotationParamH);
		printf("       L: %g\n", prms::rotationParamL);
		printf("       MC: %g\n", prms::rotationParamMc);
	}
	if (prms::usingHardeningBase)
	{
		printf(" Using basic hardening\n");
		printf("       A: %g\n", prms::hardeningParamBaseA);
		printf("       Delta: %g\n", prms::hardeningParamBaseDelta);
		printf("       Psi: %g\n", prms::hardeningParamBasePsi);
	}
	if (prms::usingHardeningBound)
	{
		printf(" Using boundary hardening\n");
		printf("       K: %g\n", prms::hardeningParamBoundK);
	}
	if (prms::debug_period > 0)
	{
		printf(" ============= DEBUG MODE =============\n");
		printf("       Period: %d\n", prms::debug_period);
		printf("       Start: %d\n", prms::DEBUG_START);
		printf("       Stop: %d\n", prms::DEBUG_STOP);
	}
	if (prms::fixedOrientations == 1)
	{
		printf(" Saving current orientations and normals\n");
	}
	if (prms::fixedOrientations == 2)
	{
		printf(" Reading saved orientations and normals\n");
	}

	Polycrystall PC;				//Создание и инициализация поликристалла
	PC.Init(prms::fragm_count);

	unsigned long t1, t2;			//Отсечки времени

	printf(" Initializing all fragments... ");
	t1 = clock();

	PC.D = prms::gradV.getSymmetryPart();
	PC.W = prms::gradV.getAntiSymmetryPart();

	std::srand(time(NULL));

	switch (prms::SurroundsGrade)			//Степень учёта соседних элементов
	{
	case 0:
	{
		prms::surroundCount = 6;			//Обычный уровень
		break;
	}
	case 1:
	{
		prms::surroundCount = 18;		//Повышенный уровень
		break;
	}
	case 2:
	{
		prms::surroundCount = 26;		//Самый высокий уровень
		break;
	}
	}

	PC.setParams();					//Заполнение всех параметров поликристалла
	PC.MakeStruct();				//Формирование фрагментной структуры
	if (prms::FRAGMENTATION) PC.MakeGrains2();//Формирование зёренной структуры


	if (prms::fixedOrientations == 2)	//Считывание записанных ориентаций
	{
		std::ifstream StreamO("DBG\\o.txt", std::ios_base::in);
		std::ifstream StreamNorm("DBG\\Norm.txt", std::ios_base::in);
		for (int q1 = 0; q1 < prms::fragm_count; q1++)
		{
			for (int q2 = 0; q2 < prms::fragm_count; q2++)
			{
				for (int q3 = 0; q3 < prms::fragm_count; q3++)
				{
					for (int i = 0; i < DIM; i++)	//Считываем значения ориентационных тензоров
					{
						for (int j = 0; j < DIM; j++)
						{
							StreamO >> PC.C[q1][q2][q3].o.C[i][j];
						}
					}
					for (int h = 0; h < prms::surroundCount; h++)//Считываем значения нормалей
					{
						for (int i = 0; i < DIM; i++)
						{
							StreamNorm >> PC.C[q1][q2][q3].normals[h].C[i];
						}
					}
				}
			}
		}
		StreamO.close();
		StreamNorm.close();
	}

	if (prms::fixedOrientations == 1)//Запоминание начальных ориентаций
	{
		std::ofstream StreamO("DBG\\o.txt", std::ios_base::out | std::ios_base::trunc);
		std::ofstream StreamNorm("DBG\\Norm.txt", std::ios_base::out | std::ios_base::trunc);
		for (int q1 = 0; q1 < prms::fragm_count; q1++)
		{
			for (int q2 = 0; q2 < prms::fragm_count; q2++)
			{
				for (int q3 = 0; q3 < prms::fragm_count; q3++)
				{
					WriteDebugInfo(StreamO, PC.C[q1][q2][q3].o.C);//Записываем значения тензоров ориентации
					for (int h = 0; h < prms::surroundCount; h++)//Записываем значения нормалей
					{
						for (int i = 0; i < DIM; i++)
						{
							StreamNorm << PC.C[q1][q2][q3].normals[h].C[i] << " ";
						}
						StreamNorm << std::endl;
					}
				}
			}
		}
		StreamO.close();
		StreamNorm.close();
	}

	if (prms::read_init_stress)	//Считывание остаточных напряжений
	{
		std::ifstream StreamSgm("DBG\\sgm.txt", std::ios_base::in);

		for (int q1 = 0; q1 < prms::fragm_count; q1++)
		{
			for (int q2 = 0; q2 < prms::fragm_count; q2++)
			{
				for (int q3 = 0; q3 < prms::fragm_count; q3++)
				{
					for (int i = 0; i < DIM; i++)	//Считываем значения тензоров
					{
						for (int j = 0; j < DIM; j++)
						{
							StreamSgm >> PC.C[q1][q2][q3].sgm.C[i][j];
						}
					}

				}
			}
		}
		StreamSgm.close();

	}
	t2 = clock();
	printf("%g sec\n", (t2 - t1) / 1000.0);

	PC.OpenFiles();				//Открытие и очистка файлов для вывода
								//Сохранение начальных полюсных фигур и ССТ
	if (prms::polus_period > 0)
	{
		printf(" Saving pole figures... ");
		t1 = clock();
		PC.SavePoleFig();
		t2 = clock();
		printf("%g sec\n", (t2 - t1) / 1000.0);
	}

	t1 = clock();		//Начальная отсечка времени
	PC.Deformate();		//Деформирование
	t2 = clock();		//Финальная отсечка времени

	if (prms::read_init_stress)	//Сохранение остаточных напряжений
	{
		PC.dbgstream[3].open("DBG\\sgm.txt", std::ios_base::out | std::ios_base::trunc);
		for (int q1 = 0; q1 < prms::fragm_count; q1++)
		{
			for (int q2 = 0; q2 < prms::fragm_count; q2++)
			{
				for (int q3 = 0; q3 < prms::fragm_count; q3++)
				{

					WriteDebugInfo(PC.dbgstream[3], PC.C[q1][q2][q3].sgm.C);

				}
			}
		}
		PC.dbgstream[3].close();
	}

	/***********************************************************
	*********       Информация о времени и шагах		********
	*********	(автоматически отправляется в панель	********
	*********	управления и программа закрывается,		********
	*********		если она была запущена из неё		********
	***********************************************************/

	if (!isnormal(PC.Strain)) printf("\n Calculation ERROR!\n");
	else printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b Done    \n");
	printf(" ==================================================\n");
	printf(" Processing time: %g sec\n", (t2 - t1) / 1000.0);
	printf(" Number of steps: %d\n", PC.CURR_STEP);
	if (!isnormal(PC.Strain))//Если не зафиксированы ошибки - закрытие
	{
		printf(" ==================================================\n");
		//printf(" Press any key or STOP button to exit...");
		system("title Done");//Меняет заголовок окна
		system("pause");
	}

	PC.CloseFiles();				//Сохранение всех полученных данных

	/************************************************************
	*******      Передача данных в панель управления      *******
	************************************************************/
	HWND hwnd;
	hwnd = ::FindWindow(NULL, (LPCTSTR)"Панель управления моделью");
	if (hwnd != NULL)
	{
		::SendMessage(hwnd, WM_USER + 1, PC.CURR_STEP, (t2 - t1));
		//Аргумент 1 - количество шагов
		//Аргумент 2 - затраченное время (мс)
	}

	return 0;
}

