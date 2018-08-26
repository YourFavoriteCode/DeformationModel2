// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"

#include <fstream>
#include <ctime>
#include <Windows.h>

#include "Params.h"
#include "Functions.h"
#include "Polycrystall.h"
#include "Loading.h"

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
	std::string paramFileName = argv[1];
	printf(" Parameters file: %s\n", paramFileName.c_str());
	if (prms::ReadParams(paramFileName.c_str()) == 1) printf(" Error in file!\n");		//Считали параметры из файла
	printf(" ==========================================\n");
	const int grainCountTotal = (int)pow(prms::grainCountLinear, 3);	//Общее кол-во фрагментов
	printf(" Fragments count: %d\n", grainCountTotal);
	printf(" Max. strain: %g\n", prms::maxStrainIntencity);
	if (prms::loadCycleCount != 1)
	{
		printf(" Number of cycles: %d\n", prms::loadCycleCount);
	}
	printf(" Integration step: %g\n", prms::dt);
	if (prms::usingRotationsHardening)
	{
		printf(" Using rotation hardening model\n");
		printf("       K1: %g\n", prms::rotationParamHardK1);
		printf("       K2: %g\n", prms::rotationParamHardK2);
	}
	if (prms::usingRotationsTaylor)
	{
		printf(" Using Taylor's rotation model\n");
	}
	if (prms::usingRotationsTrusov)
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
	if (prms::saveVariablesPeriodStep > 0)
	{
		printf(" ============= DEBUG MODE =============\n");
		printf("       Period: %d\n", prms::saveVariablesPeriodStep);
		printf("       Start: %d\n", prms::saveVariablesStartStep);
		printf("       Stop: %d\n", prms::saveVariablesStopStep);
	}
	if (prms::fixedOrientations == 1)
	{
		printf(" Saving current orientations and normals\n");
	}
	if (prms::fixedOrientations == 2)
	{
		printf(" Reading saved orientations and normals\n");
	}

	Polycrystall polycrystall;				//Создание и инициализация поликристалла
	polycrystall.init(prms::grainCountLinear);

	unsigned long t1, t2;			//Отсечки времени

	printf(" Initializing all fragments... ");
	t1 = clock();



	std::srand(time(NULL));

	switch (prms::grainSurroundGrade)			//Степень учёта соседних элементов
	{
	case prms::GRADE_BASE:
	{
		prms::grainSurroundCount = 6;			//Обычный уровень
		break;
	}
	case prms::GRADE_DETAILED:
	{
		prms::grainSurroundCount = 18;		//Повышенный уровень
		break;
	}
	case prms::GRADE_MOST_DETAILED:
	{
		prms::grainSurroundCount = 26;		//Самый высокий уровень
		break;
	}
	}

	polycrystall.setParams();					//Заполнение всех параметров поликристалла
	polycrystall.makeGrainStruct(STRUCTURE_CUBIC);				//Формирование фрагментной структуры

	if (prms::fixedOrientations == 2)	//Считывание записанных ориентаций
	{
		std::ifstream StreamO("DBG\\o.txt", std::ios_base::in);
		std::ifstream StreamNorm("DBG\\Norm.txt", std::ios_base::in);
		for (int q = 0; q < grainCountTotal; q++)
		{

			for (int i = 0; i < DIM; i++)	//Считываем значения ориентационных тензоров
			{
				for (int j = 0; j < DIM; j++)
				{
					StreamO >> polycrystall.c[q].o.c[i][j];
				}
			}
			for (int h = 0; h < prms::grainSurroundCount; h++)//Считываем значения нормалей
			{
				for (int i = 0; i < DIM; i++)
				{
					StreamNorm >> polycrystall.c[q].normals[h].c[i];
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
		for (int q = 0; q < grainCountTotal; q++)
		{
			writeDebugInfo(StreamO, polycrystall.c[q].o.c);//Записываем значения тензоров ориентации
			for (int h = 0; h < prms::grainSurroundCount; h++)//Записываем значения нормалей
			{
				for (int i = 0; i < DIM; i++)
				{
					StreamNorm << polycrystall.c[q].normals[h].c[i] << " ";
				}
				StreamNorm << std::endl;

			}
		}
		StreamO.close();
		StreamNorm.close();
	}

	if (prms::usingInititalStress)	//Считывание остаточных напряжений
	{
		std::ifstream StreamSgm("DBG\\sgm.txt", std::ios_base::in);

		for (int q = 0; q < grainCountTotal; q++)
		{
			for (int i = 0; i < DIM; i++)	//Считываем значения тензоров
			{
				for (int j = 0; j < DIM; j++)
				{
					StreamSgm >> polycrystall.c[q].sgm.c[i][j];
				}
			}
		}
		StreamSgm.close();

	}
	t2 = clock();
	printf("%g sec\n", (t2 - t1) / 1000.0);

	polycrystall.openAllFiles();				//Открытие и очистка файлов для вывода
								//Сохранение начальных полюсных фигур и ССТ
	if (prms::periodSavePolus > 0)
	{
		printf(" Saving pole figures... ");
		t1 = clock();
		polycrystall.savePoleFigData();
		t2 = clock();
		printf("%g sec\n", (t2 - t1) / 1000.0);
	}
	Loading loading;						//Нагружение
	loading.setLoad(prms::gradV, prms::loadCycleCount);

	t1 = clock();							//Начальная отсечка времени
	polycrystall.deformate(&loading);		//Деформирование
	t2 = clock();							//Финальная отсечка времени

	if (prms::usingInititalStress)			//Сохранение остаточных напряжений
	{
		polycrystall.streamDebug[3].open("DBG\\sgm.txt", std::ios_base::out | std::ios_base::trunc);
		for (int q = 0; q < grainCountTotal; q++)
		{
			writeDebugInfo(polycrystall.streamDebug[3], polycrystall.c[q].sgm.c);
		}
		polycrystall.streamDebug[3].close();
	}

	/***********************************************************
	*********       Информация о времени и шагах		********
	*********	(автоматически отправляется в панель	********
	*********	управления и программа закрывается,		********
	*********		если она была запущена из неё		********
	***********************************************************/

	if (!isNormalDouble(polycrystall.strain)) printf("\n Calculation ERROR!\n");
	else printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b Done    \n");
	printf(" ==================================================\n");
	printf(" Processing time: %g sec\n", (t2 - t1) / 1000.0);
	printf(" Number of steps: %d\n", loading.currentStep);
	if (!isNormalDouble(polycrystall.strain))//Если не зафиксированы ошибки - закрытие
	{
		printf(" ==================================================\n");
		//printf(" Press any key or STOP button to exit...");
		system("title Done");//Меняет заголовок окна
		system("pause");
	}

	polycrystall.closeAllFiles();				//Сохранение всех полученных данных

	/************************************************************
	*******      Передача данных в панель управления      *******
	************************************************************/
	HWND hwnd;
	hwnd = ::FindWindow(NULL, (LPCTSTR)"Панель управления моделью");
	if (hwnd != NULL)
	{
		::SendMessage(hwnd, WM_USER + 1, loading.currentStep, (t2 - t1));
		//Аргумент 1 - количество шагов
		//Аргумент 2 - затраченное время (мс)
	}

	return 0;
}

