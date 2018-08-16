﻿// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"
#include <cmath>

#include "Fragment.h"
#include "Params.h"
#include "Materials.h"
#include "CoreFuctions.h"

namespace model
{

	Fragment::Fragment()
	{
		//Обнуление всех переменных при инициализации
		//кроме тех, которые имеют свой конструктор
		stress = 0;
		strain = 0;
		SS_count = 0;
		size = 0;
		rot_speed = 0;
		isRotate = false;
		sum_angle = 0;
		rot_Mc = 0;
		rot_A = 0;
		rot_H = 0;
		rot_L = 0;
		rot_energy = 0;
		iter = 0;
	}

	Fragment::~Fragment()
	{

	}

	int inline H(double a)		//Функция Хэвисайда
	{
		return (a >= 0 ? 1 : 0);
	}

	void Fragment::Orientate(double a, double g,
		double cb)
	{
		/*
		* Функция задаёт ориентацию решётки на основе двух
		* углов и двух случайных чисел, генерируя равномерное
		* распределение косинуса третьего угла
		*/
		double sb = sqrt(1.0 - cb * cb);

		double sa = sin(a);
		double sg = sin(g);
		double ca = cos(a);
		double cg = cos(g);

		o.C[0][0] = ca * cg - sa * cb * sg;
		o.C[0][1] = -ca * sg - sa * cb * cg;
		o.C[0][2] = sa * sb;

		o.C[1][0] = sa * cg + ca * cb * sg;
		o.C[1][1] = -sa * sg + ca * cb * cg;
		o.C[1][2] = -ca * sb;

		o.C[2][0] = sb * sg;
		o.C[2][1] = sb * cg;
		o.C[2][2] = cb;

	}

	void Fragment::OrientateAxis(double cf, Vector a)
	{
		/*
		* Поворот решётки фрагмента вокруг
		* Заданной оси на заданный угол
		*/
		double sf = sqrt(1.0 - cf * cf);
		double sfx = sf*a.C[0];
		double sfy = sf*a.C[1];
		double sfz = sf*a.C[2];
		double COS = (1.0 - cf);

		o.C[0][0] = cf + COS * a.C[0] * a.C[0];
		o.C[0][1] = COS * a.C[0] * a.C[1] - sfz;
		o.C[0][2] = COS * a.C[0] * a.C[2] + sfy;

		o.C[1][0] = COS * a.C[0] * a.C[1] + sfz;
		o.C[1][1] = cf + COS * a.C[1] * a.C[1];
		o.C[1][2] = COS * a.C[1] * a.C[2] - sfx;

		o.C[2][0] = COS * a.C[0] * a.C[2] - sfy;
		o.C[2][1] = COS * a.C[1] * a.C[2] + sfx;
		o.C[2][2] = cf + COS * a.C[2] * a.C[2];

	}

	void Fragment::OrientateQuater(double w, double x, double y, double z)
	{
		/*
		* Ориентация решётки на основе кватерниона
		*/
		
		o.C[0][0] = 1 - 2 * y*y - 2 * z*z;
		o.C[0][1] = 2 * x*y - 2 * z*w;
		o.C[0][2] = 2 * x*z + 2 * y*w;

		o.C[1][0] = 2 * x*y + 2 * z*w;
		o.C[1][1] = 1 - 2 * x*x - 2 * z*z;
		o.C[1][2] = 2 * y*z - 2 * x*w;

		o.C[2][0] = 2 * x*z - 2 * y*w;
		o.C[2][1] = 2 * y*z + 2 * x*w;
		o.C[2][2] = 1 - 2 * x*x - 2 * y*y;

	}

	void Fragment::setMaterialParams(int material_type)
	{
		/*
		* Задание всех материальных параметров фрагмента
		* в зависимости от выбранного материала
		*/
		this->material = material_type;			//Тип материала
		int LATTICE_TYPE;						//Тип решётки
		double P1, P2, P3, P4, P5, P6;			//Упругие константы
		double c;								//Отношение c/a (для ГПУ)

		switch (material_type)//Идентификация материальных параметров
		{
		case 0://Сталь 45
		{
			SS_count = SS_COUNT_BCC;
			LATTICE_TYPE = 0;
			SS = new SlipSystem[SS_count];
			for (int i = 0; i < 24; i++)
			{
				SS[i].tc = STEEL_TC1;
			}
			for (int i = 24; i < 48; i++)
			{
				SS[i].tc = STEEL_TC2;
			}
			for (int i = 48; i < 96; i++)
			{
				SS[i].tc = STEEL_TC3;
			}
			P1 = STEEL_P1;
			P2 = STEEL_P2;
			P3 = STEEL_P3;

			break;
		}
		case 1://Медь
		{
			SS_count = SS_COUNT_FCC;
			LATTICE_TYPE = 1;
			SS = new SlipSystem[SS_count];
			for (int i = 0; i < SS_count; i++)
			{
				SS[i].tc = CUPR_TC;
			}
			P1 = CUPR_P1;
			P2 = CUPR_P2;
			P3 = CUPR_P3;

			break;
		}
		case 2://Титан
		{
			SS_count = SS_COUNT_HCP;
			LATTICE_TYPE = 2;
			SS = new SlipSystem[SS_count];
			for (int i = 0; i < 6; i++)
			{
				SS[i].tc = TITAN_TC1;
			}
			for (int i = 6; i < 12; i++)
			{
				SS[i].tc = TITAN_TC2;
			}
			for (int i = 12; i < 36; i++)
			{
				SS[i].tc = TITAN_TC3;
			}

			P1 = TITAN_P1;
			P2 = TITAN_P2;
			P3 = TITAN_P3;
			P4 = TITAN_P4;
			P5 = TITAN_P5;
			P6 = TITAN_P6;

			c = TITAN_C;

			break;
		}
		case 3://Аллюминий
		{
			SS_count = SS_COUNT_FCC;
			LATTICE_TYPE = 1;
			SS = new SlipSystem[SS_count];
			for (int i = 0; i < SS_count; i++)
			{
				SS[i].tc = 8e6;//???????
			}
			P1 = ALLUM_P1;
			P2 = ALLUM_P2;
			P3 = ALLUM_P3;

			break;
		}
		}

		if (LATTICE_TYPE == 0 || LATTICE_TYPE == 1)//Идентификация упругх констант
		{
			//Подобная симметрия допустима только для кубических решёток
			p.C[0][0][0][0] = p.C[1][1][1][1] = p.C[2][2][2][2] = P1;

			p.C[0][0][1][1] = p.C[1][1][0][0] = p.C[2][2][1][1] =
				p.C[1][1][2][2] = p.C[2][2][0][0] = p.C[0][0][2][2] = P2;

			p.C[0][1][1][0] = p.C[1][2][1][2] = p.C[0][1][0][1] =
				p.C[2][0][2][0] = p.C[0][2][0][2] = p.C[2][1][2][1] =
				p.C[1][0][1][0] = p.C[1][0][0][1] = p.C[2][1][1][2] =
				p.C[0][2][2][0] = p.C[1][2][2][1] = p.C[2][0][0][2] = P3;
		}
		else if (LATTICE_TYPE == 2)//ГПУ, альфа-титан
		{
			p.C[2][2][2][2] = P1;

			p.C[0][0][0][0] = p.C[1][1][1][1] = P2;

			p.C[0][0][2][2] = p.C[2][2][0][0] = p.C[1][1][2][2] = p.C[2][2][1][1] = P3;

			p.C[0][0][1][1] = p.C[1][1][0][0] = P4;

			p.C[1][2][1][2] = p.C[2][1][1][2] = p.C[1][2][2][1] = p.C[2][1][2][1] =
				p.C[2][0][2][0] = p.C[2][0][0][2] = p.C[0][2][0][2] = p.C[0][2][2][0] = P5;

			p.C[0][1][1][0] = p.C[1][0][1][0] = p.C[1][0][0][1] = p.C[0][1][0][1] = P6;

		}

		switch (LATTICE_TYPE)//Идентификация систем скольжения
		{
		case 0://ОЦК
		{
			//-----------[110]
			SS[0].Initialize(0, 1, 1, 1, -1, 1);
			SS[1].Initialize(0, 1, 1, 1, 1, -1);
			SS[2].Initialize(0, -1, 1, 1, 1, 1);
			SS[3].Initialize(0, -1, 1, -1, 1, 1);

			SS[4].Initialize(-1, 1, 0, 1, 1, 1);
			SS[5].Initialize(-1, 1, 0, 1, 1, -1);
			SS[6].Initialize(1, 1, 0, -1, 1, 1);
			SS[7].Initialize(1, 1, 0, 1, -1, 1);

			SS[8].Initialize(1, 0, 1, 1, 1, -1);
			SS[9].Initialize(1, 0, 1, -1, 1, 1);
			SS[10].Initialize(-1, 0, 1, 1, 1, 1);
			SS[11].Initialize(-1, 0, 1, 1, -1, 1);

			SS[12].Initialize(0, 1, 1, -1, 1, -1);
			SS[13].Initialize(0, 1, 1, -1, -1, 1);
			SS[14].Initialize(0, -1, 1, -1, -1, -1);
			SS[15].Initialize(0, -1, 1, 1, -1, -1);

			SS[16].Initialize(-1, 1, 0, -1, -1, -1);
			SS[17].Initialize(-1, 1, 0, -1, -1, 1);
			SS[18].Initialize(1, 1, 0, 1, -1, -1);
			SS[19].Initialize(1, 1, 0, -1, 1, -1);

			SS[20].Initialize(1, 0, 1, -1, -1, 1);
			SS[21].Initialize(1, 0, 1, 1, -1, -1);
			SS[22].Initialize(-1, 0, 1, -1, -1, -1);
			SS[23].Initialize(-1, 0, 1, -1, 1, -1);
			//-----------[112]
			SS[24].Initialize(2, 1, 1, 1, -1, -1);
			SS[25].Initialize(2, -1, 1, 1, 1, -1);
			SS[26].Initialize(2, -1, -1, 1, 1, 1);
			SS[27].Initialize(2, 1, -1, 1, -1, 1);

			SS[28].Initialize(1, 2, 1, -1, 1, -1);
			SS[29].Initialize(-1, 2, 1, 1, 1, -1);
			SS[30].Initialize(-1, 2, -1, 1, 1, 1);
			SS[31].Initialize(1, 2, -1, -1, 1, 1);

			SS[32].Initialize(1, 1, 2, -1, -1, 1);
			SS[33].Initialize(-1, 1, 2, 1, -1, 1);
			SS[34].Initialize(-1, -1, 2, 1, 1, 1);
			SS[35].Initialize(1, -1, 2, -1, 1, 1);

			SS[36].Initialize(2, 1, 1, -1, 1, 1);
			SS[37].Initialize(2, -1, 1, -1, -1, 1);
			SS[38].Initialize(2, -1, -1, -1, -1, -1);
			SS[39].Initialize(2, 1, -1, -1, 1, -1);

			SS[40].Initialize(1, 2, 1, 1, -1, 1);
			SS[41].Initialize(-1, 2, 1, -1, -1, 1);
			SS[42].Initialize(-1, 2, -1, -1, -1, -1);
			SS[43].Initialize(1, 2, -1, 1, -1, -1);

			SS[44].Initialize(1, 1, 2, 1, 1, -1);
			SS[45].Initialize(-1, 1, 2, -1, 1, -1);
			SS[46].Initialize(-1, -1, 2, -1, -1, -1);
			SS[47].Initialize(1, -1, 2, 1, -1, -1);
			//-----------[123]
			SS[48].Initialize(1, 2, 3, -1, -1, 1);
			SS[49].Initialize(-1, 2, 3, 1, -1, 1);
			SS[50].Initialize(-1, -2, 3, 1, 1, 1);
			SS[51].Initialize(1, -2, 3, -1, 1, 1);

			SS[52].Initialize(3, 1, 2, 1, -1, -1);
			SS[53].Initialize(3, -1, 2, 1, 1, -1);
			SS[54].Initialize(3, -1, -2, 1, 1, 1);
			SS[55].Initialize(3, 1, -2, 1, -1, 1);

			SS[56].Initialize(2, 3, 1, -1, 1, -1);
			SS[57].Initialize(2, 3, -1, -1, 1, 1);
			SS[58].Initialize(-2, 3, -1, 1, 1, 1);
			SS[59].Initialize(-2, 3, 1, 1, 1, -1);

			SS[60].Initialize(2, 1, 3, -1, -1, 1);
			SS[61].Initialize(-2, 1, 3, 1, -1, 1);
			SS[62].Initialize(-2, -1, 3, 1, 1, 1);
			SS[63].Initialize(2, -1, 3, -1, 1, 1);

			SS[64].Initialize(1, 3, 2, -1, 1, -1);
			SS[65].Initialize(-1, 3, 2, 1, 1, -1);
			SS[66].Initialize(-1, 3, -2, 1, 1, 1);
			SS[67].Initialize(1, 3, -2, -1, 1, 1);

			SS[68].Initialize(3, 2, 1, 1, -1, -1);
			SS[69].Initialize(3, -2, 1, 1, 1, -1);
			SS[70].Initialize(3, -2, -1, 1, 1, 1);
			SS[71].Initialize(3, 2, -1, 1, -1, 1);

			SS[72].Initialize(1, 2, 3, 1, 1, -1);
			SS[73].Initialize(-1, 2, 3, -1, 1, -1);
			SS[74].Initialize(-1, -2, 3, -1, -1, -1);
			SS[75].Initialize(1, -2, 3, 1, -1, -1);

			SS[76].Initialize(3, 1, 2, -1, 1, 1);
			SS[77].Initialize(3, -1, 2, -1, -1, 1);
			SS[78].Initialize(3, -1, -2, -1, -1, -1);
			SS[79].Initialize(3, 1, -2, -1, 1, -1);

			SS[80].Initialize(2, 3, 1, 1, -1, 1);
			SS[81].Initialize(2, 3, -1, 1, -1, -1);
			SS[82].Initialize(-2, 3, -1, -1, -1, -1);
			SS[83].Initialize(-2, 3, 1, -1, -1, 1);

			SS[84].Initialize(2, 1, 3, 1, 1, -1);
			SS[85].Initialize(-2, 1, 3, -1, 1, -1);
			SS[86].Initialize(-2, -1, 3, -1, -1, -1);
			SS[87].Initialize(2, -1, 3, 1, -1, -1);

			SS[88].Initialize(1, 3, 2, 1, -1, 1);
			SS[89].Initialize(-1, 3, 2, -1, -1, 1);
			SS[90].Initialize(-1, 3, -2, -1, -1, -1);
			SS[91].Initialize(1, 3, -2, 1, -1, -1);

			SS[92].Initialize(3, 2, 1, -1, 1, 1);
			SS[93].Initialize(3, -2, 1, -1, -1, 1);
			SS[94].Initialize(3, -2, -1, -1, -1, -1);
			SS[95].Initialize(3, 2, -1, -1, 1, -1);

			break;
		}
		case 1://ГЦК
		{
			SS[0].Initialize(1, -1, 1, 0, 1, 1);
			SS[1].Initialize(1, 1, -1, 0, 1, 1);
			SS[2].Initialize(1, 1, 1, 0, -1, 1);
			SS[3].Initialize(-1, 1, 1, 0, -1, 1);

			SS[4].Initialize(1, 1, 1, -1, 1, 0);
			SS[5].Initialize(1, 1, -1, -1, 1, 0);
			SS[6].Initialize(-1, 1, 1, 1, 1, 0);
			SS[7].Initialize(1, -1, 1, 1, 1, 0);

			SS[8].Initialize(1, 1, -1, 1, 0, 1);
			SS[9].Initialize(-1, 1, 1, 1, 0, 1);
			SS[10].Initialize(1, 1, 1, -1, 0, 1);
			SS[11].Initialize(1, -1, 1, -1, 0, 1);

			SS[12].Initialize(1, -1, 1, 0, -1, -1);
			SS[13].Initialize(1, 1, -1, 0, -1, -1);
			SS[14].Initialize(1, 1, 1, 0, 1, -1);
			SS[15].Initialize(-1, 1, 1, 0, 1, -1);

			SS[16].Initialize(1, 1, 1, 1, -1, 0);
			SS[17].Initialize(1, 1, -1, 1, -1, 0);
			SS[18].Initialize(-1, 1, 1, -1, -1, 0);
			SS[19].Initialize(1, -1, 1, -1, -1, 0);

			SS[20].Initialize(1, 1, -1, -1, 0, -1);
			SS[21].Initialize(-1, 1, 1, -1, 0, -1);
			SS[22].Initialize(1, 1, 1, 1, 0, -1);
			SS[23].Initialize(1, -1, 1, 1, 0, -1);

			break;
		}
		case 2://ГПУ при a/c=1.587
		{
			//Базисное скольжение
			SS[0].Initialize(0, 0, 0, 1, 1, 1, -2, 0, c);// {0001} <11-20> 
			SS[1].Initialize(0, 0, 0, 1, -1, -1, 2, 0, c);
			SS[2].Initialize(0, 0, 0, 1, -2, 1, 1, 0, c);
			SS[3].Initialize(0, 0, 0, 1, 2, -1, -1, 0, c);
			SS[4].Initialize(0, 0, 0, 1, 1, -2, 1, 0, c);
			SS[5].Initialize(0, 0, 0, 1, -1, 2, -1, 0, c);
			//Призматическое скольжение
			SS[6].Initialize(-1, 1, 0, 0, 1, 1, -2, 0, c);// {10-10} <11-20>
			SS[7].Initialize(1, -1, 0, 0, -1, -1, 2, 0, c);
			SS[8].Initialize(0, 1, -1, 0, -2, 1, 1, 0, c);
			SS[9].Initialize(0, -1, 1, 0, 2, -1, -1, 0, c);
			SS[10].Initialize(1, 0, -1, 0, 1, -2, 1, 0, c);
			SS[11].Initialize(-1, 0, 1, 0, -1, 2, -1, 0, c);
			//Пирамидальное скольжение <a + c>
			SS[12].Initialize(-1, 0, 1, 1, 1, 1, -2, 3, c);// {10-11} <11-23>
			SS[13].Initialize(1, 0, -1, 1, -1, -1, 2, 3, c);
			SS[14].Initialize(1, -1, 0, 1, -2, 1, 1, 3, c);
			SS[15].Initialize(-1, 1, 0, 1, 2, -1, -1, 3, c);
			SS[16].Initialize(0, 1, -1, 1, 1, -2, 1, 3, c);
			SS[17].Initialize(0, -1, 1, 1, -1, 2, -1, 3, c);

			SS[18].Initialize(1, 0, -1, -1, 1, 1, -2, 3, c);// {10-11} <11-23>
			SS[19].Initialize(-1, 0, 1, -1, -1, -1, 2, 3, c);
			SS[20].Initialize(-1, 1, 0, -1, -2, 1, 1, 3, c);
			SS[21].Initialize(1, -1, 0, -1, 2, -1, -1, 3, c);
			SS[22].Initialize(0, -1, 1, -1, 1, -2, 1, 3, c);
			SS[23].Initialize(0, 1, -1, -1, -1, 2, -1, 3, c);

			SS[24].Initialize(1, 0, -1, 1, 1, 1, -2, -3, c);// {10-11} <11-23>
			SS[25].Initialize(-1, 0, 1, 1, -1, -1, 2, -3, c);
			SS[26].Initialize(-1, 1, 0, 1, -2, 1, 1, -3, c);
			SS[27].Initialize(1, -1, 0, 1, 2, -1, -1, -3, c);
			SS[28].Initialize(0, -1, 1, 1, 1, -2, 1, -3, c);
			SS[29].Initialize(0, 1, -1, 1, -1, 2, -1, -3, c);

			SS[30].Initialize(-1, 0, 1, -1, 1, 1, -2, -3, c);// {10-11} <11-23>
			SS[31].Initialize(1, 0, -1, -1, -1, -1, 2, -3, c);
			SS[32].Initialize(1, -1, 0, -1, -2, 1, 1, -3, c);
			SS[33].Initialize(-1, 1, 0, -1, 2, -1, -1, -3, c);
			SS[34].Initialize(0, 1, -1, -1, 1, -2, 1, -3, c);
			SS[35].Initialize(0, -1, 1, -1, -1, 2, -1, -3, c);

			break;
		}
		}
	}

	void Fragment::NDScalc()
	{
		/*
		* Функция вычисляет все характеристики НДС фрагмента
		*/
		for (int k = 0; k < SS_count; k++)
		{
			/************************************************************
			***********              Закон Шмида             ************
			************************************************************/
			SS[k].t = 0;
			if (!prms::isSymmetrycal)
			{
				for (int i = 0; i < DIM; i++)
				{
					for (int j = 0; j < DIM; j++)
					{
						SS[k].t += sgm.C[i][j] * SS[k].n.C[i] * SS[k].b.C[j];
					}
				}
			}
			else
			{

				for (int i = 0; i < DIM; i++)
				{
					for (int j = 0; j < DIM; j++)
					{
						SS[k].t += sgm.C[i][j] * (SS[k].n.C[i] * SS[k].b.C[j] + SS[k].n.C[j] * SS[k].b.C[i]);
					}
				}
				SS[k].t /= 2.0;
			}

			/************************************************************
			***********      Соотношение Хатчинсона          ************
			************************************************************/
			/*	if (iter % 100 == 0)
				{
				FILE* G_File = fopen("Plot\\TBS.txt", "a+");

				//fprintf(G_File, "%f  \n", (SS[k].t - fabs(SS[k].tbs) - SS[k].tc));

				fclose(G_File);
				}*/
			//	SS[k].dgm = prms::dgm0 * pow((SS[k].t -(SS[k].tbs)) / SS[k].tc, prms::m) * H(SS[k].t - (SS[k].tbs) /*- SS[k].tc*/);
			SS[k].dgm = prms::dgm0 * pow(SS[k].t / SS[k].tc, prms::m) * H(SS[k].t - SS[k].tc);

			SS[k].gmm += SS[k].dgm * prms::dt; //Сдвиг по СС
		}
		/************************************************************
		**********    Вычисление неупругих деформаций     ***********
		************************************************************/
		d_in.setZero();
		if (!prms::isSymmetrycal)
		{
			for (int i = 0; i < DIM; i++)
			{
				for (int j = 0; j < DIM; j++)
				{
					for (int k = 0; k < SS_count; k++)
					{
						d_in.C[i][j] += SS[k].dgm * SS[k].n.C[i] * SS[k].b.C[j];
					}
				}
			}
		}
		else
		{
			for (int i = 0; i < DIM; i++)
			{
				for (int j = 0; j < DIM; j++)
				{
					for (int k = 0; k < SS_count; k++)

					{
						d_in.C[i][j] += SS[k].dgm * (SS[k].n.C[i] * SS[k].b.C[j] + SS[k].n.C[j] * SS[k].b.C[i]);
					}
				}
			}
			d_in /= 2.0;
		}


		/************************************************************
		**********               Закон Гука               ***********
		************************************************************/
		//Коротационные производные
		Tensor R = om*sgm - sgm*om;						//Связанная со спином
		//Tensor J = w*sgm-sgm*w;						//Яуманна-Нолла
		//Tensor CR = (w + d)*sgm + sgm*Transp(w + d);	//Коттер-Ривлина
		//Tensor OL = Transp(w + d)*sgm + sgm*(w + d);	//Олдройда

		dsgm.setZero();
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int l = 0; l < DIM; l++)
				{
					for (int n = 0; n < DIM; n++)
					{
						dsgm.C[i][j] += p.C[i][j][l][n] * (d.C[n][l] - d_in.C[n][l]);
					}
				}
			}
		}
		dsgm += R;
		//dsgm += J;
		//dsgm -= OL;
		//dsgm += CR;

		/************************************************************
		**********  Интенсивности деформаций и напряжений  **********
		************************************************************/

		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				e.C[i][j] += d.C[i][j] * prms::dt;
				sgm.C[i][j] += dsgm.C[i][j] * prms::dt;
			}
		}

		strain = e.doubleScalMult(e);
		stress = sgm.doubleScalMult(sgm);

		strain = SQRT2_3 * sqrt(strain);
		stress = SQRT3_2 * sqrt(stress);
	}

	double Fragment::DisorientMeasure(int h)
	{
		double M = 0;
		for (int i = 0; i < 3; i++)
		{
			Vector e;	//Направления [100], [110], [111] в КСК
			for (int j = 0; j <= i; j++) e.C[j] = 1;
			e.Normalize();

			Vector e1 = ScalMult(e, o);//Перевод вектора в ЛСК
			Vector e2 = ScalMult(e, surrounds[h].o);
			e1.Normalize();
			e2.Normalize();
			//Теперь нужно отразить оба эти вектора в одну четверть сферы (x>0,y>0,z>0)
			for (int k = 0; k < 3; k++)
			{
				e1.C[k] = fabs(e1.C[k]);
				e2.C[k] = fabs(e2.C[k]);
			}
			//А потом свернуть ещё пополам
			if (e1.C[0] < e1.C[1])
			{
				double buf = e1.C[0];
				e1.C[0] = e1.C[1];
				e1.C[1] = buf;
			}
			if (e2.C[0] < e2.C[1])
			{
				double buf = e2.C[0];
				e2.C[0] = e2.C[1];
				e2.C[1] = buf;
			}

			double teta1 = atan(sqrt(e1.C[0] * e1.C[0] + e1.C[1] * e1.C[1]) / e1.C[2]);
			double fi1 = atan(e1.C[1] / e1.C[0]);

			double teta2 = atan(sqrt(e2.C[0] * e2.C[0] + e2.C[1] * e2.C[1]) / e2.C[2]);
			double fi2 = atan(e2.C[1] / e2.C[0]);
			//Нахождение длины дуги между двумя точками на единичной сфере
			double L = acos(cos(teta1)*cos(teta2) + sin(teta1)*sin(teta2)*cos(fi1 - fi2)) / PI_2;//Нормировка по PI/2
			M = (L > M) ? L : M;
		}
		return M;
	}
}