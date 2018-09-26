// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"

#include <cmath>
#include <ctime>
#include <omp.h>
#include <fstream>
#include "vector"

#include "Polycrystall.h"
#include "Params.h"
#include "Distributions.h"
#include "Tension.h"
#include "Rotations.h"
#include "Hardening.h"
#include "Functions.h"
#include "GrainStructure.h"
#include "Fragmentation.h"

namespace model
{
	Polycrystall::Polycrystall()
	{
		strain = 0;
		stress = 0;
		stepShowProgress = 0;
		stepSavePoleFig = 0;
		stepSavePlot = 0;
		stepSaveInternalData = 0;
		periodShowProgress = 400;
		fileCount = 16;
	}

	Polycrystall::~Polycrystall()
	{

	}

	void Polycrystall::openAllFiles()
	{
		truncPoleFigFiles();				//Очистка всех файлов полюсных фигур
		truncSSTFiles();

		streamDataX = new std::ofstream[20];//Открытие файлов для записи кривых НДС
		streamDataY = new std::ofstream[20];
		streamInternalVars = new std::ofstream[1];

		//Файлы для вывода макро-данных
		if (prms::saveMacro && prms::saveIntensity)
		{
			streamDataX[0].open("Plot\\macroXint.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[0].open("Plot\\macroYint.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}

		if (prms::saveMacro && prms::save11)
		{
			streamDataX[1].open("Plot\\macroX11.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[1].open("Plot\\macroY11.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMacro && prms::save12)
		{
			streamDataX[2].open("Plot\\macroX12.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[2].open("Plot\\macroY12.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMacro && prms::save13)
		{
			streamDataX[3].open("Plot\\macroX13.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[3].open("Plot\\macroY13.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMacro && prms::save21)
		{
			streamDataX[4].open("Plot\\macroX21.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[4].open("Plot\\macroY21.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMacro && prms::save22)
		{
			streamDataX[5].open("Plot\\macroX22.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[5].open("Plot\\macroY22.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMacro && prms::save23)
		{
			streamDataX[6].open("Plot\\macroX23.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[6].open("Plot\\macroY23.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMacro && prms::save31)
		{
			streamDataX[7].open("Plot\\macroX31.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[7].open("Plot\\macroY31.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMacro && prms::save32)
		{
			streamDataX[8].open("Plot\\macroX32.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[8].open("Plot\\macroY32.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMacro && prms::save33)
		{
			streamDataX[9].open("Plot\\macroX33.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[9].open("Plot\\macroY33.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}

		//Файлы для мезо-данных

		if (prms::saveMeso && prms::saveIntensity)
		{
			streamDataX[10].open("Plot\\mesoXint.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[10].open("Plot\\mesoYint.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}

		if (prms::saveMeso && prms::save11)
		{
			streamDataX[11].open("Plot\\mesoX11.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[11].open("Plot\\mesoY11.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMeso && prms::save12)
		{
			streamDataX[12].open("Plot\\mesoX12.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[12].open("Plot\\mesoY12.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMeso && prms::save13)
		{
			streamDataX[13].open("Plot\\mesoX13.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[13].open("Plot\\mesoY13.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMeso && prms::save21)
		{
			streamDataX[14].open("Plot\\mesoX21.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[14].open("Plot\\mesoY21.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMeso && prms::save22)
		{
			streamDataX[15].open("Plot\\mesoX22.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[15].open("Plot\\mesoY22.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMeso && prms::save23)
		{
			streamDataX[16].open("Plot\\mesoX23.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[16].open("Plot\\mesoY23.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMeso && prms::save31)
		{
			streamDataX[17].open("Plot\\mesoX31.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[17].open("Plot\\mesoY31.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMeso && prms::save32)
		{
			streamDataX[18].open("Plot\\mesoX32.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[18].open("Plot\\mesoY32.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}
		if (prms::saveMeso && prms::save33)
		{
			streamDataX[19].open("Plot\\mesoX33.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
			streamDataY[19].open("Plot\\mesoY33.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
		}

		streamInternalVars[0].open("Plot\\ActiveSS.dat", std::ios_base::out | std::ios_base::trunc | std::ios::binary);

		streamDebug = new std::ofstream[fileCount];
		if (prms::saveVariablesPeriodStep > 0)				//Открытие файлов для отладочных данных
		{
			streamDebug[0].open("DBG\\o.txt", std::ios_base::out | std::ios_base::trunc);
			streamDebug[1].open("DBG\\e.txt", std::ios_base::out | std::ios_base::trunc);
			streamDebug[2].open("DBG\\d.txt", std::ios_base::out | std::ios_base::trunc);
			streamDebug[3].open("DBG\\sgm.txt", std::ios_base::out | std::ios_base::trunc);
			streamDebug[4].open("DBG\\om.txt", std::ios_base::out | std::ios_base::trunc);
			streamDebug[5].open("DBG\\dsgm.txt", std::ios_base::out | std::ios_base::trunc);
			streamDebug[6].open("DBG\\din.txt", std::ios_base::out | std::ios_base::trunc);
			streamDebug[7].open("DBG\\w.txt", std::ios_base::out | std::ios_base::trunc);
			streamDebug[8].open("DBG\\dgamma.txt", std::ios_base::out | std::ios_base::trunc);
			streamDebug[9].open("DBG\\t.txt", std::ios_base::out | std::ios_base::trunc);
			streamDebug[10].open("DBG\\Macro_D.txt", std::ios_base::out | std::ios_base::trunc);
			streamDebug[11].open("DBG\\Macro_Din.txt", std::ios_base::out | std::ios_base::trunc);
			streamDebug[12].open("DBG\\Macro_Sgm.txt", std::ios_base::out | std::ios_base::trunc);
			streamDebug[13].open("DBG\\Macro_dSgm.txt", std::ios_base::out | std::ios_base::trunc);
			streamDebug[14].open("DBG\\Macro_E.txt", std::ios_base::out | std::ios_base::trunc);
			streamDebug[15].open("DBG\\VOL_M.txt", std::ios_base::out | std::ios_base::trunc);
		}

		streamDataTest = new std::ofstream[6];
		streamDataTest[0].open("Test0.txt", std::ios_base::out | std::ios_base::trunc);
		streamDataTest[1].open("Test1.txt", std::ios_base::out | std::ios_base::trunc);
		streamDataTest[2].open("Test2.txt", std::ios_base::out | std::ios_base::trunc);
		streamDataTest[3].open("Test3.txt", std::ios_base::out | std::ios_base::trunc);
		streamDataTest[4].open("Test4.txt", std::ios_base::out | std::ios_base::trunc);
		streamDataTest[5].open("Test5.txt", std::ios_base::out | std::ios_base::trunc);
	}

	void Polycrystall::closeAllFiles()
	{
		for (int i = 0; i < 20; i++)
		{
			streamDataX[i].close();
			streamDataY[i].close();
		}
		streamInternalVars[0].close();
		for (int i = 0; i < 6; i++)
		{
			streamDataTest[i].close();
		}

		if (prms::saveVariablesPeriodStep > 0)
		{
			for (int i = 0; i < fileCount; i++)
			{
				streamDebug[i].close();
			}
		}
	}

	void Polycrystall::init(int count)
	{
		grainCount = count;
		totalGrainCount = (int)pow(count, 3);
		c = std::vector<Grain>(totalGrainCount); //Выделение памяти под массив
	}

	Grain* Polycrystall::findGrainById(int id)
	{
		for (Grain g : c)
		{
			if (g.position == id)
			{
				return &g;
			}
		}
	}

	int Polycrystall::findGrainPosById(int id)
	{
		for (int i = 0; i < c.size(); i++)
		{
			if (c[i].position == id)
			{
				return i;
			}
		}
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
				c[q].setMaterialParams(prms::materialType);
			}
			else
			{
				c[q].setMaterialParams(another_material);
			}

			c[q].rotationParamMc = prms::rotationParamMc;	//Раздача начальных критических моментов
			c[q].rotationParamA = prms::rotationParamA;	//и параметров модели ротаций
			c[q].rotationParamH = prms::rotationParamH;
			c[q].rotationParamL = prms::rotationParamL;
			c[q].position = q;//Получение порядкового номера фрагмента

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
					c[q].orientate(a, g, cb);
				}
				else//Получение ориентационного тензора (КСК=ЛСК)
				{
					c[q].o.setUnit();
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
					c[q].orientateAxis(cf, axis);
				}
				else//Получение ориентационного тензора (КСК=ЛСК)
				{
					c[q].o.setUnit();
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
					c[q].orientateQuater(w, x, y, z);
				}
				else//Получение ориентационного тензора (КСК=ЛСК)
				{
					c[q].o.setUnit();
				}
			}

			//Задание размеров фрагментов
			switch (prms::grainSizeDistribLaw)
			{
			case prms::DISTRIB_UNIFORM:
			{
				c[q].size = UniformDistrib(prms::grainSizeDistribM, prms::grainSizeDistribD);
				break;
			}
			case prms::DISTRIB_NORMAL:
			{
				c[q].size = NormalDistrib(prms::grainSizeDistribM, prms::grainSizeDistribD);
				break;
			}
			case prms::DISTRIB_LOGNORMAL:
			{
				c[q].size = LogNormalDistrib(prms::grainSizeDistribM, prms::grainSizeDistribD);
				break;
			}
			case prms::DISTRIB_EXPONENT:
			{
				c[q].size = ExpDistrib(prms::grainSizeDistribM);//Только один параметр
				break;
			}
			}			

		}
	}
	
	void Polycrystall::makeGrainStruct() {
	
		structure = new GrainStructure(this, true);
		structure->makeVoronoiStructure();

	}

	//TODO: Перенос работы с файлами в отдельный класс

	void Polycrystall::savePoleFigData()
	{
		for (int q = 0; q < totalGrainCount; q++)
		{
			getPoleFig(&c[q]);
			if (prms::usingStandardTriangleSaving) getSST(&c[q]);

		}
	}

	void Polycrystall::saveDebugData()
	{
		for (int i = 0; i < fileCount; i++)//Визуальное разделение шагов 
		{
			streamDebug[i] << "#########################      STEP " << loading->currentStep << "      #########################" << std::endl << std::endl;
		}

		//Запись тензоров каждого из зерен или фрагментов
		for (int q = 0; q < totalGrainCount; q++)
		{

			writeDebugInfo(streamDebug[0], c[q].o.c);
			writeDebugInfo(streamDebug[1], c[q].e.c);
			writeDebugInfo(streamDebug[2], c[q].d.c);
			writeDebugInfo(streamDebug[3], c[q].sgm.c);
			writeDebugInfo(streamDebug[4], c[q].om.c);
			writeDebugInfo(streamDebug[5], c[q].dsgm.c);
			writeDebugInfo(streamDebug[6], c[q].d_in.c);
			writeDebugInfo(streamDebug[7], c[q].w.c);
			for (int f = 0; f < c[q].ssCount; f++)
			{
				streamDebug[8] << c[q].ss[f].dgm << " ";
			}
			streamDebug[8] << std::endl << std::endl;
			for (int f = 0; f < c[q].ssCount; f++)
			{
				streamDebug[9] << c[q].ss[f].t << " ";
			}
			streamDebug[9] << std::endl << std::endl;
			streamDebug[15] << c[q].rotationMoment.c[0] << " " << c[q].rotationMoment.c[1] << " " << c[q].rotationMoment.c[2] << std::endl;

		}
		//Запись тензоров представительного объема
		writeDebugInfo(streamDebug[10], D.c);
		writeDebugInfo(streamDebug[11], D_in.c);
		writeDebugInfo(streamDebug[12], Sgm.c);
		writeDebugInfo(streamDebug[13], dSgm.c);
		writeDebugInfo(streamDebug[14], E.c);
	}

	void Polycrystall::load(bool unload)
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
				D_in += c[q].d_in;
				P += c[q].p.toLsk(c[q].o);
			}

			D_in /= (totalGrainCount);
			P /= (totalGrainCount);

			//Симметризация тензора упругих констант
			P.symmetrize();

			D = !unload ? TensionStrainCalc(P, D_in, D.c[0][0]) : UnloadingStrainCalc(P, D_in, Sgm, loading->paramLambda);

			strain = SQRT2_3 * sqrt(E.doubleScalMult(E));//Вычисление интенсивности деформаций

			dSgm = TensionStressCalc(P, D_in, D);
			//dSgm *= prms::dt;				//Приращение напряжений на шаге
			Sgm += dSgm * prms::dt;
			stress = SQRT3_2 * sqrt(Sgm.doubleScalMult(Sgm));//Вычисление интенсивности напряжений
		}
		else
		{
			stress = 0;		//Вычисление интенсивностей осреднением
			strain = 0;
			for (int q = 0; q < totalGrainCount; q++)
			{
				strain += c[q].strain;
				stress += c[q].stress;

			}
			strain /= totalGrainCount;
			stress /= totalGrainCount;
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

			Tensor O = c[q].o;
			Tensor OT = O;
			OT.transp();
			c[q].d = O * D*OT;//Гипотеза Фойгта
			c[q].w = O * W*OT /*- C[q].om*/;//Расширенная

			c[q].sgm = O * c[q].sgm*OT;
			c[q].d_in = O * c[q].d_in*OT;


			/***************************************************
			***********       Пересчитываем НДС      ***********
			***************************************************/

			c[q].stressStrainCalc();

			if (prms::usingHardeningBase)			//Базовое упрочнение
			{
				hardeningBase(&c[q]);
			}

			if (prms::usingRotationsTaylor)		//Ротации по Тейлору
			{
				rotateByTaylor(&c[q]);
			}

			if (prms::usingRotationsTrusov && prms::usingRotationsHardening)	//Ротационное упрочнение
			{
				rotationHardening(&c[q]);
			}

			if (prms::usingHardeningBound)	//Зернограничное упрочнение
			{
				hardeningBoundary(&c[q]);
			}

			if (prms::usingRotationsTrusov)		//Ротации по Трусову
			{
				rotateByTrusov(&c[q]);
			}


			/**************************************************
			************       Переходим в ЛСК       **********
			**************************************************/

			c[q].sgm = OT * c[q].sgm*O;
			c[q].d_in = OT * c[q].d_in*O;

		}

		if (prms::usingFragmentation)		// Фрагментация зерен
		{
			structure->fragmentate();
		}
		
		if (!prms::trueUniaxial && !unload)		//Этот блок нужен исключительно для работы с энергией!
		{
			Sgm.setZero();
			D_in.setZero();
			for (int q = 0; q < totalGrainCount; q++)
			{

				Sgm += c[q].sgm;
				D_in += c[q].d_in;

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
			progress = prms::trueUniaxial ? fabs(E.c[0][0]) : strain;
			progress = progress / prms::maxStrainIntencity * 100.0;

			if (!(loading->cycleCount == 1 || loading->cycle == 0))	//Для многоцикловых нагружений
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
			progress = loading->maxStress / fabs(Sgm.c[0][0]) * 100.0;	//Индикация прогресса при разгрузке
		}

		int period = unload ? periodShowProgress / 40 : periodShowProgress;

		if (stepShowProgress == period)
		{
			stepShowProgress = 0;
			//Курсор двигается на 6 символов влево, записывается новое значение с точностью 5.2 и знак %
			printf("\b\b\b\b\b\b%0*.*f%%", 5, 2, progress);
		}

		/************************************************************
		***********	    Запись данных для графиков НДС    ***********
		************************************************************/

		if ((progress - stepSavePlot > prms::periodSavePlot || unload) && prms::periodSavePlot > 0)
		{

			if (prms::saveMacro)	//Запись компонент тензоров макроуровня
			{
				if (prms::saveIntensity)
				{
					streamDataX[0].write((char *)&strain, sizeof(double));
					streamDataY[0].write((char *)&stress, sizeof(double));
				}
				if (prms::save11)
				{
					streamDataX[1].write((char *)&E.c[0][0], sizeof(double));
					streamDataY[1].write((char *)&Sgm.c[0][0], sizeof(double));
				}
				if (prms::save12)
				{
					streamDataX[2].write((char *)&E.c[0][1], sizeof(double));
					streamDataY[2].write((char *)&Sgm.c[0][1], sizeof(double));
				}
				if (prms::save13)
				{
					streamDataX[3].write((char *)&E.c[0][2], sizeof(double));
					streamDataY[3].write((char *)&Sgm.c[0][2], sizeof(double));
				}
				if (prms::save21)
				{
					streamDataX[4].write((char *)&E.c[1][0], sizeof(double));
					streamDataY[4].write((char *)&Sgm.c[1][0], sizeof(double));
				}
				if (prms::save22)
				{
					streamDataX[5].write((char *)&E.c[1][1], sizeof(double));
					streamDataY[5].write((char *)&Sgm.c[1][1], sizeof(double));
				}
				if (prms::save23)
				{
					streamDataX[6].write((char *)&E.c[1][2], sizeof(double));
					streamDataY[6].write((char *)&Sgm.c[1][2], sizeof(double));
				}
				if (prms::save31)
				{
					streamDataX[7].write((char *)&E.c[2][0], sizeof(double));
					streamDataY[7].write((char *)&Sgm.c[2][0], sizeof(double));
				}
				if (prms::save32)
				{
					streamDataX[8].write((char *)&E.c[2][1], sizeof(double));
					streamDataY[8].write((char *)&Sgm.c[2][1], sizeof(double));
				}
				if (prms::save33)
				{
					streamDataX[9].write((char *)&E.c[2][2], sizeof(double));
					streamDataY[9].write((char *)&Sgm.c[2][2], sizeof(double));
				}
			}

			if (prms::saveMeso)	//Запись компонент тензоров мезоуровня
			{
				for (int q = 0; q < totalGrainCount; q++)
				{

					if (prms::saveIntensity)
					{
						streamDataX[10].write((char *)&c[q].strain, sizeof(double));
						streamDataY[10].write((char *)&c[q].stress, sizeof(double));
					}
					if (prms::save11)
					{
						streamDataX[11].write((char *)&c[q].e.c[0][0], sizeof(double));
						streamDataY[11].write((char *)&c[q].sgm.c[0][0], sizeof(double));
					}
					if (prms::save12)
					{
						streamDataX[12].write((char *)&c[q].e.c[0][1], sizeof(double));
						streamDataY[12].write((char *)&c[q].sgm.c[0][1], sizeof(double));
					}
					if (prms::save13)
					{
						streamDataX[13].write((char *)&c[q].e.c[0][2], sizeof(double));
						streamDataY[13].write((char *)&c[q].sgm.c[0][2], sizeof(double));
					}
					if (prms::save21)
					{
						streamDataX[14].write((char *)&c[q].e.c[1][0], sizeof(double));
						streamDataY[14].write((char *)&c[q].sgm.c[1][0], sizeof(double));
					}
					if (prms::save22)
					{
						streamDataX[15].write((char *)&c[q].e.c[1][1], sizeof(double));
						streamDataY[15].write((char *)&c[q].sgm.c[1][1], sizeof(double));
					}
					if (prms::save23)
					{
						streamDataX[16].write((char *)&c[q].e.c[1][2], sizeof(double));
						streamDataY[16].write((char *)&c[q].sgm.c[1][2], sizeof(double));
					}
					if (prms::save31)
					{
						streamDataX[17].write((char *)&c[q].e.c[2][0], sizeof(double));
						streamDataY[17].write((char *)&c[q].sgm.c[2][0], sizeof(double));
					}
					if (prms::save32)
					{
						streamDataX[18].write((char *)&c[q].e.c[2][1], sizeof(double));
						streamDataY[18].write((char *)&c[q].sgm.c[2][1], sizeof(double));
					}
					if (prms::save33)
					{
						streamDataX[19].write((char *)&c[q].e.c[2][2], sizeof(double));
						streamDataY[19].write((char *)&c[q].sgm.c[2][2], sizeof(double));
					}

				}
			}

			double ActiveSysCount = 0;			//Среднее кол-во активных систем скольжения на шаге
			double RotEnergy = 0;				//Энергия ротаций на шаге
			double RotSpeed = 0;				//Средняя скорость вращения на шаге
			int RotCount = 0;					//Кол-во вращающихся фрагментов
			double Mc = 0;
			double angle = 0;
			double H = 0;
			Vector M;
			for (int q = 0; q < totalGrainCount; q++)
			{
				for (int i = 0; i < c[q].ssCount; i++)
				{
					if (c[q].ss[i].dgm > EPS) ActiveSysCount++;//Подсчёт активных СС
				}

				if (c[q].isRotate) RotCount++;		//Подсчёт вращающихся решёток
				RotEnergy += c[q].rotationEnergy;		//Суммирование энергий вращения
				RotSpeed += c[q].rotationSpeed;		//Суммирование скоростей вращения
				Mc += c[q].rotationParamMc;
				angle += c[q].rotationTotalAngle;
				H += c[q].rotationParamH;
				M += c[q].rotationMoment;
			}
			M /= totalGrainCount;
			H /= totalGrainCount;
			angle /= totalGrainCount;
			Mc /= totalGrainCount;
			ActiveSysCount /= totalGrainCount;
			if (prms::saveActiveSS) streamInternalVars[0].write((char *)&ActiveSysCount, sizeof ActiveSysCount);//Запись кол-ва активных СС
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
			streamDataTest[0] << RotCount << std::endl;
			streamDataTest[1] << RotSpeed << std::endl;
			streamDataTest[2] << RotEnergy << std::endl;
			streamDataTest[3] << StepEnergy << std::endl;
			streamDataTest[4] << StepEnergy_in << std::endl;
			streamDataTest[5] << M.getNorm() << std::endl;
			stepSavePlot = progress;
		}

		/************************************************************
		***********	      Сохранение полюсных фигур	      ***********
		************************************************************/
		if (progress - stepSavePoleFig > prms::periodSavePolus && prms::periodSavePolus > 0)
		{
			savePoleFigData();
			stepSavePoleFig = progress;
		}

		/************************************************************
		***********	       Запись пошаговых данных	      ***********
		************************************************************/
		if (loading->currentStep >= prms::saveVariablesStartStep && loading->currentStep <= prms::saveVariablesStopStep && stepSaveInternalData == prms::saveVariablesPeriodStep)
		{
			stepSaveInternalData = 0;
			saveDebugData();
		}
		loading->currentStep++;
		stepShowProgress++;
		if (loading->currentStep >= prms::saveVariablesStartStep && loading->currentStep <= prms::saveVariablesStopStep)
		{
			stepSaveInternalData++;
		}
	}

	void Polycrystall::deformate(Loading *newLoading)
	{
		/****************************************
		Функция деформирования представительного
		объема поликристала
		*****************************************/
		this->loading = newLoading;
		D = loading->gradV.getSymmetryPart();
		W = loading->gradV.getAntiSymmetryPart();

		omp_set_num_threads(prms::ompThreadCount);	//Кол-во используемых потоков

		if (prms::trueUniaxial)
		{
			/***************************************************
			*Выбор растягивающей компоненты тензора D.
			*Относительно неё на каждом шаге будет решаться СЛАУ
			***************************************************/
			loading->tensionComponent = D.c[0][0];
		}
		for (loading->cycle = 0; loading->cycle < loading->cycleCount; loading->cycle++)
		{
			unsigned long t1, t2;

			stepSavePlot = 0;		//Обнуление счетчиков для периодического вывода данных в файлы
			stepSavePoleFig = 0;
			stepShowProgress = 0;
			if (loading->cycleCount > 1)
			{
				printf("\n Loading #%d", loading->cycle + 1);
				t1 = clock();		//Начальная отсечка времени
			}
			printf("\n        00.00%");
			double* counter;
			counter = !prms::trueUniaxial ? &strain : &E.c[0][0];
			//Для циклических нагружений цикл ведется по значению растягивающей компоненты
			//а для остальных - по значению интенсивности тензора деформации

			while (fabs(*counter) < prms::maxStrainIntencity)	//Цикл по деформациям
			{
				load(false);
			}
			if (loading->cycleCount > 1)
			{
				t2 = clock();		//Конечная отсечка времени
				printf("\b\b\b\b\b\bDone in %g sec", (t2 - t1) / 1000.0);
			}
			//Если нужнен второй этап нагружения (не циклика), то его нужно добавлять здесь

			if (prms::withUnloading)	//Упругая разгрузка
			{
				printf("\n Unloading #%d", loading->cycle + 1);
				printf("\n        00.00%");
				t1 = clock();		//Начальная отсечка времени
				while (fabs(Sgm.c[0][0]) > loading->maxStress) //Цикл по напряжениям
				{
					load(true);
				}
				t2 = clock();		//Конечная отсечка времени
				printf("\b\b\b\b\b\bDone in %g sec", (t2 - t1) / 1000.0);
			}

			if (loading->cycleCount > 1)	//Циклическое знакопеременное нагружение
			{
				D.c[0][0] = pow(-1, loading->cycle + 1) * loading->tensionComponent;	//Меняем знак растягивающей компоненты
				prms::maxStrainIntencity += prms::maxStrainIntencity * loading->additionalStrain;	//Повышаем предел интенсивности
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