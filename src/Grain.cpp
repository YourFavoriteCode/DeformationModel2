// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"
#include <cmath>

#include "Grain.h"
#include "Params.h"
#include "Materials.h"

namespace model
{
	Grain::Grain()
	{
		stress = 0;
		strain = 0;
		ssCount = 0;
		size = 0;
		rotationSpeed = 0;
		isRotate = false;
		rotationTotalAngle = 0;
		rotationParamMc = 0;
		rotationParamA = 0;
		rotationParamH = 0;
		rotationParamL = 0;
		rotationEnergy = 0;
		position = -1;
	}

	Grain::~Grain()
	{

	}

	//Функция Хэвисайда
	int inline heavisideFunc(double a)		
	{
		return (a >= 0 ? 1 : 0);
	}

	/*
	* Функция задаёт ориентацию решётки на основе двух
	* углов и двух случайных чисел, генерируя равномерное
	* распределение косинуса третьего угла
	*/
	void Grain::orientate(double a, double g, double cb)
	{
		double sb = sqrt(1.0 - cb * cb);
		double sa = sin(a);
		double sg = sin(g);
		double ca = cos(a);
		double cg = cos(g);

		o.c[0][0] = ca * cg - sa * cb * sg;
		o.c[0][1] = -ca * sg - sa * cb * cg;
		o.c[0][2] = sa * sb;

		o.c[1][0] = sa * cg + ca * cb * sg;
		o.c[1][1] = -sa * sg + ca * cb * cg;
		o.c[1][2] = -ca * sb;

		o.c[2][0] = sb * sg;
		o.c[2][1] = sb * cg;
		o.c[2][2] = cb;

	}

	/*
	* Поворот решётки фрагмента вокруг
	* Заданной оси на заданный угол
	*/
	void Grain::orientateAxis(double cf, Vector a)
	{
		double sf = sqrt(1.0 - cf * cf);
		double sfx = sf*a.c[0];
		double sfy = sf*a.c[1];
		double sfz = sf*a.c[2];
		double cosf = (1.0 - cf);

		o.c[0][0] = cosf * a.c[0] * a.c[0] + cf;
		o.c[0][1] = cosf * a.c[0] * a.c[1] - sfz;
		o.c[0][2] = cosf * a.c[0] * a.c[2] + sfy;

		o.c[1][0] = cosf * a.c[1] * a.c[0] + sfz;
		o.c[1][1] = cosf * a.c[1] * a.c[1] + cf;
		o.c[1][2] = cosf * a.c[1] * a.c[2] - sfx;

		o.c[2][0] = cosf * a.c[2] * a.c[0] - sfy;
		o.c[2][1] = cosf * a.c[2] * a.c[1] + sfx;
		o.c[2][2] = cosf * a.c[2] * a.c[2] + cf;

	}

	/*
	* Ориентация решётки на основе кватерниона
	*/
	void Grain::orientateQuater(double w, double x, double y, double z)
	{
		o.c[0][0] = 1 - 2 * y*y - 2 * z*z;
		o.c[0][1] = 2 * x*y - 2 * z*w;
		o.c[0][2] = 2 * x*z + 2 * y*w;

		o.c[1][0] = 2 * x*y + 2 * z*w;
		o.c[1][1] = 1 - 2 * x*x - 2 * z*z;
		o.c[1][2] = 2 * y*z - 2 * x*w;

		o.c[2][0] = 2 * x*z - 2 * y*w;
		o.c[2][1] = 2 * y*z + 2 * x*w;
		o.c[2][2] = 1 - 2 * x*x - 2 * y*y;

	}

	/*
	* Задание всех материальных параметров фрагмента
	* в зависимости от выбранного материала
	*/
	void Grain::setMaterialParams(int materialType)
	{
		this->materialType = materialType;		// Тип материала
		int latticeType;						// Тип решётки
		double p1, p2, p3, p4, p5, p6;			// Упругие константы
		double c;								// Отношение c/a (для ГПУ)

		switch (materialType)//Идентификация материальных параметров
		{
		case 0://Сталь 45
		{
			ssCount = SS_COUNT_BCC;
			latticeType = 0;
			ss = new SlipSystem[ssCount];
			for (int i = 0; i < 24; i++)
			{
				ss[i].tc = STEEL_TC1;
			}
			for (int i = 24; i < 48; i++)
			{
				ss[i].tc = STEEL_TC2;
			}
			for (int i = 48; i < 96; i++)
			{
				ss[i].tc = STEEL_TC3;
			}
			p1 = STEEL_P1;
			p2 = STEEL_P2;
			p3 = STEEL_P3;

			break;
		}
		case 1://Медь
		{
			ssCount = SS_COUNT_FCC;
			latticeType = 1;
			ss = new SlipSystem[ssCount];
			for (int i = 0; i < ssCount; i++)
			{
				ss[i].tc = CUPR_TC;
			}
			p1 = CUPR_P1;
			p2 = CUPR_P2;
			p3 = CUPR_P3;

			break;
		}
		case 2://Титан
		{
			ssCount = SS_COUNT_HCP;
			latticeType = 2;
			ss = new SlipSystem[ssCount];
			for (int i = 0; i < 6; i++)
			{
				ss[i].tc = TITAN_TC1;
			}
			for (int i = 6; i < 12; i++)
			{
				ss[i].tc = TITAN_TC2;
			}
			for (int i = 12; i < 36; i++)
			{
				ss[i].tc = TITAN_TC3;
			}

			p1 = TITAN_P1;
			p2 = TITAN_P2;
			p3 = TITAN_P3;
			p4 = TITAN_P4;
			p5 = TITAN_P5;
			p6 = TITAN_P6;

			c = TITAN_C;

			break;
		}
		case 3://Аллюминий
		{
			ssCount = SS_COUNT_FCC;
			latticeType = 1;
			ss = new SlipSystem[ssCount];
			for (int i = 0; i < ssCount; i++)
			{
				ss[i].tc = 8e6;//???????
			}
			p1 = ALLUM_P1;
			p2 = ALLUM_P2;
			p3 = ALLUM_P3;

			break;
		}
		}

		if (latticeType == 0 || latticeType == 1)//Идентификация упругх констант
		{
			//Подобная симметрия допустима только для кубических решёток
			p.c[0][0][0][0] = p.c[1][1][1][1] = p.c[2][2][2][2] = p1;

			p.c[0][0][1][1] = p.c[1][1][0][0] = p.c[2][2][1][1] =
				p.c[1][1][2][2] = p.c[2][2][0][0] = p.c[0][0][2][2] = p2;

			p.c[0][1][1][0] = p.c[1][2][1][2] = p.c[0][1][0][1] =
				p.c[2][0][2][0] = p.c[0][2][0][2] = p.c[2][1][2][1] =
				p.c[1][0][1][0] = p.c[1][0][0][1] = p.c[2][1][1][2] =
				p.c[0][2][2][0] = p.c[1][2][2][1] = p.c[2][0][0][2] = p3;
		}
		else if (latticeType == 2)//ГПУ, альфа-титан
		{
			p.c[2][2][2][2] = p1;

			p.c[0][0][0][0] = p.c[1][1][1][1] = p2;

			p.c[0][0][2][2] = p.c[2][2][0][0] = p.c[1][1][2][2] = p.c[2][2][1][1] = p3;

			p.c[0][0][1][1] = p.c[1][1][0][0] = p4;

			p.c[1][2][1][2] = p.c[2][1][1][2] = p.c[1][2][2][1] = p.c[2][1][2][1] =
				p.c[2][0][2][0] = p.c[2][0][0][2] = p.c[0][2][0][2] = p.c[0][2][2][0] = p5;

			p.c[0][1][1][0] = p.c[1][0][1][0] = p.c[1][0][0][1] = p.c[0][1][0][1] = p6;

		}

		switch (latticeType)//Идентификация систем скольжения
		{
		case 0://ОЦК
		{
			//-----------[110]
			ss[0].Initialize(0, 1, 1, 1, -1, 1);
			ss[1].Initialize(0, 1, 1, 1, 1, -1);
			ss[2].Initialize(0, -1, 1, 1, 1, 1);
			ss[3].Initialize(0, -1, 1, -1, 1, 1);

			ss[4].Initialize(-1, 1, 0, 1, 1, 1);
			ss[5].Initialize(-1, 1, 0, 1, 1, -1);
			ss[6].Initialize(1, 1, 0, -1, 1, 1);
			ss[7].Initialize(1, 1, 0, 1, -1, 1);

			ss[8].Initialize(1, 0, 1, 1, 1, -1);
			ss[9].Initialize(1, 0, 1, -1, 1, 1);
			ss[10].Initialize(-1, 0, 1, 1, 1, 1);
			ss[11].Initialize(-1, 0, 1, 1, -1, 1);

			ss[12].Initialize(0, 1, 1, -1, 1, -1);
			ss[13].Initialize(0, 1, 1, -1, -1, 1);
			ss[14].Initialize(0, -1, 1, -1, -1, -1);
			ss[15].Initialize(0, -1, 1, 1, -1, -1);

			ss[16].Initialize(-1, 1, 0, -1, -1, -1);
			ss[17].Initialize(-1, 1, 0, -1, -1, 1);
			ss[18].Initialize(1, 1, 0, 1, -1, -1);
			ss[19].Initialize(1, 1, 0, -1, 1, -1);

			ss[20].Initialize(1, 0, 1, -1, -1, 1);
			ss[21].Initialize(1, 0, 1, 1, -1, -1);
			ss[22].Initialize(-1, 0, 1, -1, -1, -1);
			ss[23].Initialize(-1, 0, 1, -1, 1, -1);
			//-----------[112]
			ss[24].Initialize(2, 1, 1, 1, -1, -1);
			ss[25].Initialize(2, -1, 1, 1, 1, -1);
			ss[26].Initialize(2, -1, -1, 1, 1, 1);
			ss[27].Initialize(2, 1, -1, 1, -1, 1);

			ss[28].Initialize(1, 2, 1, -1, 1, -1);
			ss[29].Initialize(-1, 2, 1, 1, 1, -1);
			ss[30].Initialize(-1, 2, -1, 1, 1, 1);
			ss[31].Initialize(1, 2, -1, -1, 1, 1);

			ss[32].Initialize(1, 1, 2, -1, -1, 1);
			ss[33].Initialize(-1, 1, 2, 1, -1, 1);
			ss[34].Initialize(-1, -1, 2, 1, 1, 1);
			ss[35].Initialize(1, -1, 2, -1, 1, 1);

			ss[36].Initialize(2, 1, 1, -1, 1, 1);
			ss[37].Initialize(2, -1, 1, -1, -1, 1);
			ss[38].Initialize(2, -1, -1, -1, -1, -1);
			ss[39].Initialize(2, 1, -1, -1, 1, -1);

			ss[40].Initialize(1, 2, 1, 1, -1, 1);
			ss[41].Initialize(-1, 2, 1, -1, -1, 1);
			ss[42].Initialize(-1, 2, -1, -1, -1, -1);
			ss[43].Initialize(1, 2, -1, 1, -1, -1);

			ss[44].Initialize(1, 1, 2, 1, 1, -1);
			ss[45].Initialize(-1, 1, 2, -1, 1, -1);
			ss[46].Initialize(-1, -1, 2, -1, -1, -1);
			ss[47].Initialize(1, -1, 2, 1, -1, -1);
			//-----------[123]
			ss[48].Initialize(1, 2, 3, -1, -1, 1);
			ss[49].Initialize(-1, 2, 3, 1, -1, 1);
			ss[50].Initialize(-1, -2, 3, 1, 1, 1);
			ss[51].Initialize(1, -2, 3, -1, 1, 1);

			ss[52].Initialize(3, 1, 2, 1, -1, -1);
			ss[53].Initialize(3, -1, 2, 1, 1, -1);
			ss[54].Initialize(3, -1, -2, 1, 1, 1);
			ss[55].Initialize(3, 1, -2, 1, -1, 1);

			ss[56].Initialize(2, 3, 1, -1, 1, -1);
			ss[57].Initialize(2, 3, -1, -1, 1, 1);
			ss[58].Initialize(-2, 3, -1, 1, 1, 1);
			ss[59].Initialize(-2, 3, 1, 1, 1, -1);

			ss[60].Initialize(2, 1, 3, -1, -1, 1);
			ss[61].Initialize(-2, 1, 3, 1, -1, 1);
			ss[62].Initialize(-2, -1, 3, 1, 1, 1);
			ss[63].Initialize(2, -1, 3, -1, 1, 1);

			ss[64].Initialize(1, 3, 2, -1, 1, -1);
			ss[65].Initialize(-1, 3, 2, 1, 1, -1);
			ss[66].Initialize(-1, 3, -2, 1, 1, 1);
			ss[67].Initialize(1, 3, -2, -1, 1, 1);

			ss[68].Initialize(3, 2, 1, 1, -1, -1);
			ss[69].Initialize(3, -2, 1, 1, 1, -1);
			ss[70].Initialize(3, -2, -1, 1, 1, 1);
			ss[71].Initialize(3, 2, -1, 1, -1, 1);

			ss[72].Initialize(1, 2, 3, 1, 1, -1);
			ss[73].Initialize(-1, 2, 3, -1, 1, -1);
			ss[74].Initialize(-1, -2, 3, -1, -1, -1);
			ss[75].Initialize(1, -2, 3, 1, -1, -1);

			ss[76].Initialize(3, 1, 2, -1, 1, 1);
			ss[77].Initialize(3, -1, 2, -1, -1, 1);
			ss[78].Initialize(3, -1, -2, -1, -1, -1);
			ss[79].Initialize(3, 1, -2, -1, 1, -1);

			ss[80].Initialize(2, 3, 1, 1, -1, 1);
			ss[81].Initialize(2, 3, -1, 1, -1, -1);
			ss[82].Initialize(-2, 3, -1, -1, -1, -1);
			ss[83].Initialize(-2, 3, 1, -1, -1, 1);

			ss[84].Initialize(2, 1, 3, 1, 1, -1);
			ss[85].Initialize(-2, 1, 3, -1, 1, -1);
			ss[86].Initialize(-2, -1, 3, -1, -1, -1);
			ss[87].Initialize(2, -1, 3, 1, -1, -1);

			ss[88].Initialize(1, 3, 2, 1, -1, 1);
			ss[89].Initialize(-1, 3, 2, -1, -1, 1);
			ss[90].Initialize(-1, 3, -2, -1, -1, -1);
			ss[91].Initialize(1, 3, -2, 1, -1, -1);

			ss[92].Initialize(3, 2, 1, -1, 1, 1);
			ss[93].Initialize(3, -2, 1, -1, -1, 1);
			ss[94].Initialize(3, -2, -1, -1, -1, -1);
			ss[95].Initialize(3, 2, -1, -1, 1, -1);

			break;
		}
		case 1://ГЦК
		{
			ss[0].Initialize(1, -1, 1, 0, 1, 1);
			ss[1].Initialize(1, 1, -1, 0, 1, 1);
			ss[2].Initialize(1, 1, 1, 0, -1, 1);
			ss[3].Initialize(-1, 1, 1, 0, -1, 1);

			ss[4].Initialize(1, 1, 1, -1, 1, 0);
			ss[5].Initialize(1, 1, -1, -1, 1, 0);
			ss[6].Initialize(-1, 1, 1, 1, 1, 0);
			ss[7].Initialize(1, -1, 1, 1, 1, 0);

			ss[8].Initialize(1, 1, -1, 1, 0, 1);
			ss[9].Initialize(-1, 1, 1, 1, 0, 1);
			ss[10].Initialize(1, 1, 1, -1, 0, 1);
			ss[11].Initialize(1, -1, 1, -1, 0, 1);

			ss[12].Initialize(1, -1, 1, 0, -1, -1);
			ss[13].Initialize(1, 1, -1, 0, -1, -1);
			ss[14].Initialize(1, 1, 1, 0, 1, -1);
			ss[15].Initialize(-1, 1, 1, 0, 1, -1);

			ss[16].Initialize(1, 1, 1, 1, -1, 0);
			ss[17].Initialize(1, 1, -1, 1, -1, 0);
			ss[18].Initialize(-1, 1, 1, -1, -1, 0);
			ss[19].Initialize(1, -1, 1, -1, -1, 0);

			ss[20].Initialize(1, 1, -1, -1, 0, -1);
			ss[21].Initialize(-1, 1, 1, -1, 0, -1);
			ss[22].Initialize(1, 1, 1, 1, 0, -1);
			ss[23].Initialize(1, -1, 1, 1, 0, -1);

			break;
		}
		case 2://ГПУ при a/c=1.587
		{
			//Базисное скольжение
			// {0001} <11-20> 
			ss[0].Initialize(0, 0, 0, 1, 1, 1, -2, 0, c);
			ss[1].Initialize(0, 0, 0, 1, -1, -1, 2, 0, c);
			ss[2].Initialize(0, 0, 0, 1, -2, 1, 1, 0, c);
			ss[3].Initialize(0, 0, 0, 1, 2, -1, -1, 0, c);
			ss[4].Initialize(0, 0, 0, 1, 1, -2, 1, 0, c);
			ss[5].Initialize(0, 0, 0, 1, -1, 2, -1, 0, c);
			// Призматическое скольжение
			// {10-10} <11-20>
			ss[6].Initialize(-1, 1, 0, 0, 1, 1, -2, 0, c);
			ss[7].Initialize(1, -1, 0, 0, -1, -1, 2, 0, c);
			ss[8].Initialize(0, 1, -1, 0, -2, 1, 1, 0, c);
			ss[9].Initialize(0, -1, 1, 0, 2, -1, -1, 0, c);
			ss[10].Initialize(1, 0, -1, 0, 1, -2, 1, 0, c);
			ss[11].Initialize(-1, 0, 1, 0, -1, 2, -1, 0, c);
			// Пирамидальное скольжение <a + c>
			// {10-11} <11-23>
			ss[12].Initialize(-1, 0, 1, 1, 1, 1, -2, 3, c);
			ss[13].Initialize(1, 0, -1, 1, -1, -1, 2, 3, c);
			ss[14].Initialize(1, -1, 0, 1, -2, 1, 1, 3, c);
			ss[15].Initialize(-1, 1, 0, 1, 2, -1, -1, 3, c);
			ss[16].Initialize(0, 1, -1, 1, 1, -2, 1, 3, c);
			ss[17].Initialize(0, -1, 1, 1, -1, 2, -1, 3, c);
			// {10-11} <11-23>
			ss[18].Initialize(1, 0, -1, -1, 1, 1, -2, 3, c);
			ss[19].Initialize(-1, 0, 1, -1, -1, -1, 2, 3, c);
			ss[20].Initialize(-1, 1, 0, -1, -2, 1, 1, 3, c);
			ss[21].Initialize(1, -1, 0, -1, 2, -1, -1, 3, c);
			ss[22].Initialize(0, -1, 1, -1, 1, -2, 1, 3, c);
			ss[23].Initialize(0, 1, -1, -1, -1, 2, -1, 3, c);
			// {10-11} <11-23>
			ss[24].Initialize(1, 0, -1, 1, 1, 1, -2, -3, c);
			ss[25].Initialize(-1, 0, 1, 1, -1, -1, 2, -3, c);
			ss[26].Initialize(-1, 1, 0, 1, -2, 1, 1, -3, c);
			ss[27].Initialize(1, -1, 0, 1, 2, -1, -1, -3, c);
			ss[28].Initialize(0, -1, 1, 1, 1, -2, 1, -3, c);
			ss[29].Initialize(0, 1, -1, 1, -1, 2, -1, -3, c);
			// {10-11} <11-23>
			ss[30].Initialize(-1, 0, 1, -1, 1, 1, -2, -3, c);
			ss[31].Initialize(1, 0, -1, -1, -1, -1, 2, -3, c);
			ss[32].Initialize(1, -1, 0, -1, -2, 1, 1, -3, c);
			ss[33].Initialize(-1, 1, 0, -1, 2, -1, -1, -3, c);
			ss[34].Initialize(0, 1, -1, -1, 1, -2, 1, -3, c);
			ss[35].Initialize(0, -1, 1, -1, -1, 2, -1, -3, c);

			break;
		}
		}
	}

	/*
	* Функция вычисляет все характеристики НДС фрагмента
	*/
	void Grain::stressStrainCalc()
	{
		d_in.setZero();

		for (int k = 0; k < ssCount; k++)
		{
			/************************************************************
			***********              Закон Шмида             ************
			************************************************************/
			
			// Ориентационный тензор системы скольжения
			Tensor o = ss[k].o;
			if (prms::isSymmetrycal)
			{
				o = o.getSymmetryPart();
			}
			ss[k].t = sgm.doubleScalMult(o);

			/************************************************************
			***********      Соотношение Хатчинсона          ************
			************************************************************/
			
			//	ss[k].dgm = prms::shearRateLawDgm0 * pow((ss[k].t -(ss[k].tbs)) / ss[k].tc, prms::m) * heavisideFunc(ss[k].t - (ss[k].tbs) /*- ss[k].tc*/);
			ss[k].dgm = prms::shearRateLawDgm0 * pow(ss[k].t / ss[k].tc, prms::shearRateLawM) * heavisideFunc(ss[k].t - ss[k].tc);

			ss[k].gmm += ss[k].dgm * prms::dt;
			/************************************************************
			**********    Вычисление неупругих деформаций     ***********
			************************************************************/
			d_in += ss[k].dgm * o;
		}

		/************************************************************
		**********               Закон Гука               ***********
		************************************************************/
		dsgm = p.doubleScalMult(d - d_in);

		// Коротационные производные
		dsgm += om*sgm - sgm*om;					// Связанная со спином
		//dsgm += w*sgm-sgm*w;						// Яуманна-Нолла
		//dsgm += (w + d)*sgm + sgm*Transp(w + d);	// Коттер-Ривлина
		//dsgm -= Transp(w + d)*sgm + sgm*(w + d);	// Олдройда

		/************************************************************
		**********  Интенсивности деформаций и напряжений  **********
		************************************************************/
		e += d * prms::dt;
		sgm += dsgm * prms::dt;

		strain = e.doubleScalMult(e);
		stress = sgm.doubleScalMult(sgm);

		strain = SQRT2_3 * sqrt(strain);
		stress = SQRT3_2 * sqrt(stress);
	}

	double Grain::disorientMeasure(int h)
	{
		double M = 0;
		for (int i = 0; i < 3; i++)
		{
			Vector e;	//Направления [100], [110], [111] в КСК
			for (int j = 0; j <= i; j++) e.c[j] = 1;
			e.normalize();

			Vector e1 = scalMult(e, o);//Перевод вектора в ЛСК
			Vector e2 = scalMult(e, neighbors[h]->o);
			e1.normalize();
			e2.normalize();
			//Теперь нужно отразить оба эти вектора в одну четверть сферы (x>0,y>0,z>0)
			for (int k = 0; k < 3; k++)
			{
				e1.c[k] = fabs(e1.c[k]);
				e2.c[k] = fabs(e2.c[k]);
			}
			//А потом свернуть ещё пополам
			if (e1.c[0] < e1.c[1])
			{
				double buf = e1.c[0];
				e1.c[0] = e1.c[1];
				e1.c[1] = buf;
			}
			if (e2.c[0] < e2.c[1])
			{
				double buf = e2.c[0];
				e2.c[0] = e2.c[1];
				e2.c[1] = buf;
			}

			double teta1 = atan(sqrt(e1.c[0] * e1.c[0] + e1.c[1] * e1.c[1]) / e1.c[2]);
			double fi1 = atan(e1.c[1] / e1.c[0]);

			double teta2 = atan(sqrt(e2.c[0] * e2.c[0] + e2.c[1] * e2.c[1]) / e2.c[2]);
			double fi2 = atan(e2.c[1] / e2.c[0]);
			//Нахождение длины дуги между двумя точками на единичной сфере
			double L = acos(cos(teta1)*cos(teta2) + sin(teta1)*sin(teta2)*cos(fi1 - fi2)) / PI_2;//Нормировка по PI/2
			M = (L > M) ? L : M;
		}
		return M;
	}

	void Grain::clearTopology()
	{
		neighbors.clear();
		areas.clear();
		normals.clear();
	}
}