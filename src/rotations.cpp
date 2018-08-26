﻿// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"
#include <fstream>

#include "Rotations.h"
#include "Params.h"
#include "Fragment.h"
#include "Distributions.h"

namespace model
{
	void Rotate(Fragment* f, double dFi, const Vector a)
	{
		/*
		* Поворот решётки фрагмента вокруг
		* Заданной оси на заданный угол
		*/
		double CosFi = cos(dFi);
		double SinFiX = sin(dFi)*a.c[0];
		double SinFiY = sin(dFi)*a.c[1];
		double SinFiZ = sin(dFi)*a.c[2];
		double COS = (1.0 - CosFi);

		Tensor dO;			//Тензор поворота на шаге

		dO.c[0][0] = CosFi + COS * a.c[0] * a.c[0];
		dO.c[0][1] = COS * a.c[0] * a.c[1] - SinFiZ;
		dO.c[0][2] = COS * a.c[0] * a.c[2] + SinFiY;

		dO.c[1][0] = COS * a.c[0] * a.c[1] + SinFiZ;
		dO.c[1][1] = CosFi + COS * a.c[1] * a.c[1];
		dO.c[1][2] = COS * a.c[1] * a.c[2] - SinFiX;

		dO.c[2][0] = COS * a.c[0] * a.c[2] - SinFiY;
		dO.c[2][1] = COS * a.c[1] * a.c[2] + SinFiX;
		dO.c[2][2] = CosFi + COS * a.c[2] * a.c[2];

		Tensor buf = dO*f->o;
		f->o = buf;
	}
	
	void Taylor_rotations(Fragment *f)
	{
		f->om.setZero();
		//f->om = f->w - f->d_in.getAntiSymmetryPart();
		for (int i = 0; i < DIM; i++)//Спин решётки
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < f->SS_count; k++)
				{
					f->om.c[i][j] -= f->SS[k].dgm * (f->SS[k].n.c[i] * f->SS[k].b.c[j] - f->SS[k].b.c[i] * f->SS[k].n.c[j]);
				}
			}
		}
		f->om /= 2.0;
		f->om += f->w;
		
		Vector e;				//Ось вращения решётки
		e.set(f->om.c[1][2], f->om.c[2][0], f->om.c[0][1]);
		
		
		double dFi = e.getNorm();
		f->isRotate = dFi > EPS*1000;
		if (f->isRotate)
		{
			f->rot_speed = dFi;
			dFi *= prms::dt;			//Получаем угол поворота
			f->sum_angle += dFi;
			e.normalize();
			Rotate(f, dFi, e);
	
		}
		else
		{
			f->rot_speed = 0;
		
		}
	}

	void Rotation_hardening(Fragment *f)
	{
		if (f->sum_angle > 100*EPS)
		{
			double dmc = prms::rotationParamHardK1 + prms::rotationParamHardK2*f->sum_angle / f->volume;
			f->rot_Mc += dmc*prms::dt;//Приращение критического момента
		}
	}

	void Trusov_rotations(Fragment *f)
	{
		Vector dM;						//Производная вектор-момента
		double S = f->size*f->size;		//Площадь фасетки (полная)
		for (int h = 0; h < prms::grainSurroundCount; h++)//Пробегаем по всем соседям фрагмента
		{
			if (f->contact[h] == 0) continue;//Если нет контакта - пропускаем
			Tensor d_in1, d_in2;
			for (int i = 0; i < DIM; i++)
			{
				for (int j = 0; j < DIM; j++)
				{
					for (int k = 0; k < f->SS_count; k++)
					{
						if (f->SS[k].b.scalMult(f->normals[h]) < 0) continue; //Скольжение от границы - вклад не вносится
						d_in1.c[i][j] += f->SS[k].dgm * (f->SS[k].n.c[i] * f->SS[k].b.c[j]);
					}
					for (int k = 0; k < f->neighbors[h]->SS_count; k++)
					{
						//if (f->SS[k].b.ScalMult(f->normals[h]) < 0) continue; //Скольжение от границы - вклад не вносится
						d_in2.c[i][j] += f->neighbors[h]->SS[k].dgm * (f->neighbors[h]->SS[k].n.c[i] * f->neighbors[h]->SS[k].b.c[j]);
					}
				}
			}
			Tensor Lp = d_in1 - d_in2;
			Lp.transp();
			Tensor buf = VectMult(f->normals[h], Lp);
			Vector dm = ScalMult(buf, f->normals[h]);//Поверхностный вектор-момент
			dm *= f->rot_L;
			Vector b1 = ScalMult(f->om, dm);//(коротационная производная)
			Vector b2 = ScalMult(dm, f->om);
			dm = dm + b1 - b2;
			double c;		//Определяет площадь контакта (в долях от полной площади стороны куба)
			if (h < 6) c = UniformDistrib(0.93, 0.07);
			else if (h < 14) c =  UniformDistrib(0.1, 0.05);
			else c = UniformDistrib(0.01, 0.005);
			dM += dm*S*c;
		}
		dM /= f->volume;
		double dMnorm = dM.getNorm();
		Vector M = f->moment + dM*prms::dt;
		f->moment = M;
		double norm = M.getNorm();
		if (norm > f->rot_Mc || norm == -1)
		{
			norm = f->rot_Mc;
		}
	
	
		double pr = M.scalMult(dM);
	
		double dFi = 0;		//Скорость вращения
		if (norm == f->rot_Mc && pr >= 0)	//Пластические и упругие развороты
		{
			dFi = f->rot_A * dMnorm + f->rot_H * norm;
		}
		else
		{
			dFi = f->rot_A * dMnorm;		//Только упругие развороты
		}
		
		f->isRotate = (dFi > EPS*1e4);
		if (f->isRotate)
		{
			Vector e = M;					//Ось вращения решётки сонаправлена с вектором момента
			e.normalize();
			f->rot_speed = dFi;
			dFi *= prms::dt;
			f->sum_angle += dFi;		//Накопленный угол вращения увеличивается
			Rotate(f, dFi, e);			//Вращение решётки
			f->rot_energy = norm*dFi;	//Считаем энергию ротаций

			f->om.setZero();			//Спин решётки
			for (int i = 0; i < DIM; i++)
			{
				for (int j = 0; j < DIM; j++)
				{
					for (int k = 0; k < DIM; k++)
					{
						f->om.c[i][j] -= LeviCivit(i, j, k) * e.c[k] * f->rot_speed;
					}
				}
			}
			
		}
		else
		{
			f->om.setZero();
			f->rot_speed = 0;		//Решётка не вращается
			f->rot_energy = 0;		//Энергия вращения равна нулю
		}
	}

	void inline SavePoints(Tensor O, const char *file, const int i, const int j, const int k)
	{
		Vector e;
		e.set(i, j, k);
		e.normalize();
		Vector e1 = ScalMult(e, O);//Перевод в ЛСК
		e1.normalize();
		std::ofstream of;
		of.open(file, std::ios::out | std::ios_base::app | std::ios::binary);
		//Запись в виде x,y,z
		of.write((char *)&e1.c[0], sizeof(double));
		of.write((char *)&e1.c[1], sizeof(double));
		of.write((char *)&e1.c[2], sizeof(double));
		of.close();
	}

	void GetPoleFig(Fragment *f)
	{
		/*---Семейство направлений [001]---*/
		SavePoints(f->o, "Polus\\S001.dat", 0, 0, 1);
		SavePoints(f->o, "Polus\\S010.dat", 0, 1, 0);
		SavePoints(f->o, "Polus\\S100.dat", 1, 0, 0);

		/*---Семейство направлений [011]---*/
		SavePoints(f->o, "Polus\\S011.dat", 0, 1, 1);
		SavePoints(f->o, "Polus\\S110.dat", 1, 1, 0);
		SavePoints(f->o, "Polus\\S101.dat", 1, 0, 1);

		/*---Семейство направлений [111]---*/
		SavePoints(f->o, "Polus\\S1-11.dat", 1, -1, 1);
		SavePoints(f->o, "Polus\\S-111.dat", -1, 1, 1);
		SavePoints(f->o, "Polus\\S111.dat", 1, 1, 1);
	}

	void inline GetSSTPoint(Tensor O, double dFi, const char *FileName, const int i, const int j, const int k)
	{
		Vector a;
		a.set(i, j, k);
		a.normalize();

		double CosFi = cos(dFi);
		double SinFiX = sin(dFi) * a.c[0];
		double SinFiY = sin(dFi) * a.c[1];
		double SinFiZ = sin(dFi) * a.c[2];
		double COS = (1.0 - CosFi);

		Tensor dO;			//Тензор поворота на шаге

		dO.c[0][0] = CosFi + COS * a.c[0] * a.c[0];
		dO.c[0][1] = COS * a.c[0] * a.c[1] - SinFiZ;
		dO.c[0][2] = COS * a.c[0] * a.c[2] + SinFiY;

		dO.c[1][0] = COS * a.c[0] * a.c[1] + SinFiZ;
		dO.c[1][1] = CosFi + COS * a.c[1] * a.c[1];
		dO.c[1][2] = COS * a.c[1] * a.c[2] - SinFiX;

		dO.c[2][0] = COS * a.c[0] * a.c[2] - SinFiY;
		dO.c[2][1] = COS * a.c[1] * a.c[2] + SinFiX;
		dO.c[2][2] = CosFi + COS * a.c[2] * a.c[2];
		
		Tensor R = dO * O;
		Vector v = ScalMult(a, R);
		v.normalize();
	
		std::ofstream Of;
		Of.open(FileName, std::ios::out | std::ios_base::app | std::ios::binary);
		//Запись координат в формате [xyz] подряд для всех зёрен 
		Of.write((char *)&v.c[0], sizeof(double));
		Of.write((char *)&v.c[1], sizeof(double));
		Of.write((char *)&v.c[2], sizeof(double));
		Of.close();
	}

	void GetSST(Fragment *f)
	{
		for (double fi = 0; fi < 2 * PI; fi += PI_2)//Ось четвёртого порядка
		{
			GetSSTPoint(f->o, fi, "Polus\\SST001.dat", 0, 0, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST001.dat", 0, 1, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST001.dat", 1, 0, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST001.dat", 0, 0, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST001.dat", 0, -1, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST001.dat", -1, 0, 0);
		}
		for (double fi = 0; fi < 2 * PI; fi += PI)//Ось второго порядка
		{
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 0, 1, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 0, 1, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 1, 0, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 1, 1, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 0, -1, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", -1, 1, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", -1, 0, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", -1, 0, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 1, 0, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 1, -1, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", -1, -1, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 0, -1, -1);
		}
		for (double fi = 0; fi < 2 * PI; fi += 2.0 * PI / 3.0)//Ось третьего порядка
		{
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", 1, 1, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", 1, 1, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", 1, -1, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", -1, 1, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", 1, -1, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", -1, 1, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", -1, -1, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", -1, -1, -1);
		}

	}
}