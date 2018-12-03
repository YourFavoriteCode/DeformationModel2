// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"
#include <fstream>

#include "Rotations.h"
#include "Params.h"
#include "Grain.h"
#include "Distributions.h"

namespace model
{
	void rotate(Grain* f, double dFi, const Vector a)
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

		Tensor buf = dO * f->o;
		f->o = buf;
	}

	void rotateByTaylor(Grain *f)
	{
		f->om.setZero();
		//f->om = f->w - f->d_in.getAntiSymmetryPart();
		for (int i = 0; i < DIM; i++)//Спин решётки
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < f->ssCount; k++)
				{
					f->om.c[i][j] -= f->ss[k].dgm * (f->ss[k].n.c[i] * f->ss[k].b.c[j] - f->ss[k].b.c[i] * f->ss[k].n.c[j]);
				}
			}
		}
		f->om /= 2.0;
		f->om += f->w;

		Vector e;				//Ось вращения решётки
		e.set(f->om.c[1][2], f->om.c[2][0], f->om.c[0][1]);


		double dFi = e.getNorm();
		f->isRotate = dFi > EPS * 1000;
		if (f->isRotate)
		{
			f->rotationSpeed = dFi;
			dFi *= prms::dt;			//Получаем угол поворота
			f->rotationTotalAngle += dFi;
			e.normalize();
			rotate(f, dFi, e);

		}
		else
		{
			f->rotationSpeed = 0;

		}
	}

	void rotationHardening(Grain *f)
	{
		if (f->rotationTotalAngle > 100 * EPS)
		{
			double dmc = prms::rotationParamHardK1 + prms::rotationParamHardK2*f->rotationTotalAngle / f->volume;
			f->rotationParamMc += dmc * prms::dt;//Приращение критического момента
		}
	}

	void rotateByTrusov(Grain *f)
	{
		// Производная вектор-момента
		Vector dM;						
		// Цикл по всем фасеткам данного зерна
		for (int h = 0; h < f->neighbors.size(); h++)
		{
			Tensor d_in1, d_in2;
			//for (int i = 0; i < DIM; i++)
			//{
			//	for (int j = 0; j < DIM; j++)
			//	{
			//		for (int k = 0; k < f->ssCount; k++)
			//		{
			//			if (f->ss[k].b.scalMult(f->normals[h]) < 0) continue; //Скольжение от границы - вклад не вносится
			//			d_in1.c[i][j] += f->ss[k].dgm * (f->ss[k].n.c[i] * f->ss[k].b.c[j]);
			//		}
			//		for (int k = 0; k < f->neighbors[h]->ssCount; k++)
			//		{
			//			//if (f->ss[k].b.scalMult(f->normals[h]) < 0) continue; //Скольжение от границы - вклад не вносится
			//			d_in2.c[i][j] += f->neighbors[h]->ss[k].dgm * (f->neighbors[h]->ss[k].n.c[i] * f->neighbors[h]->ss[k].b.c[j]);
			//		}
			//	}
			//}

			for (int k = 0; k < f->ssCount; k++)
			{
				d_in1 += f->ss[k].dgm * f->ss[k].o;
			}
			for (int k = 0; k < f->neighbors[h]->ssCount; k++)
			{
				d_in2 += f->neighbors[h]->ss[k].dgm * f->neighbors[h]->ss[k].o;
			}
			// Скачок тензора пластической деформации на границе зерен
			Tensor Lp = Transp(d_in1 - d_in2);
			Tensor buf = vectMult(f->normals[h], Lp);
			// Поверхностный вектор-момент, действующий на фасетку
			Vector dm = scalMult(buf, f->normals[h]);
			dm *= f->rotationParamL;
			// Слагаемые для коротационной производной
			Vector b1 = scalMult(f->om, dm);
			Vector b2 = scalMult(dm, f->om);
			dm = dm + b1 - b2;
			dM += dm * f->areas[h];
		}
		dM /= f->volume;
		double dMnorm = dM.getNorm();
		f->rotationMoment += dM * prms::dt;
		Vector M = f->rotationMoment;
		double norm = M.getNorm();
		if (norm > f->rotationParamMc || norm == -1)
		{
			norm = f->rotationParamMc;
		}

		double pr = M.scalMult(dM);
		// Вычисление скорости вращения
		// Только упругая составляющая разворотов
		double dFi = f->rotationParamA * dMnorm;
		if (norm == f->rotationParamMc && pr >= 0)
		{
			// Пластическая составляющая
			dFi += f->rotationParamH * norm;
		}

		f->isRotate = (dFi > EPS);
		if (f->isRotate)
		{
			// Ось вращения решётки сонаправлена с вектором момента
			Vector e = M;
			e.normalize();
			f->rotationSpeed = dFi;
			dFi *= prms::dt;
			// Накопленный угол вращения увеличивается
			f->rotationTotalAngle += dFi;
			// Вращение решетки
			rotate(f, dFi, e);
			// Энергия ротаций
			f->rotationEnergy = norm * dFi;
			// Спин решётки
			f->om.setZero();
			for (int i = 0; i < DIM; i++)
			{
				for (int j = 0; j < DIM; j++)
				{
					for (int k = 0; k < DIM; k++)
					{
						f->om.c[i][j] -= LeviCivit(i, j, k) * e.c[k] * f->rotationSpeed;
					}
				}
			}

		}
		else
		{
			f->om.setZero();
			f->rotationSpeed = 0;		//Решетка не вращается
			f->rotationEnergy = 0;		//Энергия вращения равна нулю
		}
	}

	void inline savePoints(Tensor O, const char *file, const int i, const int j, const int k)
	{
		Vector e;
		e.set(i, j, k);
		e.normalize();
		Vector e1 = scalMult(e, O);//Перевод в ЛСК
		e1.normalize();
		std::ofstream of;
		of.open(file, std::ios::out | std::ios_base::app | std::ios::binary);
		//Запись в виде x,y,z
		of.write((char *)&e1.c[0], sizeof(double));
		of.write((char *)&e1.c[1], sizeof(double));
		of.write((char *)&e1.c[2], sizeof(double));
		of.close();
	}

	void getPoleFig(Grain *f)
	{
		/*---Семейство направлений [001]---*/
		savePoints(f->o, "Polus\\S001.dat", 0, 0, 1);
		savePoints(f->o, "Polus\\S010.dat", 0, 1, 0);
		savePoints(f->o, "Polus\\S100.dat", 1, 0, 0);

		/*---Семейство направлений [011]---*/
		savePoints(f->o, "Polus\\S011.dat", 0, 1, 1);
		savePoints(f->o, "Polus\\S110.dat", 1, 1, 0);
		savePoints(f->o, "Polus\\S101.dat", 1, 0, 1);

		/*---Семейство направлений [111]---*/
		savePoints(f->o, "Polus\\S1-11.dat", 1, -1, 1);
		savePoints(f->o, "Polus\\S-111.dat", -1, 1, 1);
		savePoints(f->o, "Polus\\S111.dat", 1, 1, 1);
	}

	void inline saveSSTPoints(Tensor O, double dFi, const char *FileName, const int i, const int j, const int k)
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
		Vector v = scalMult(a, R);
		v.normalize();

		std::ofstream Of;
		Of.open(FileName, std::ios::out | std::ios_base::app | std::ios::binary);
		//Запись координат в формате [xyz] подряд для всех зёрен 
		Of.write((char *)&v.c[0], sizeof(double));
		Of.write((char *)&v.c[1], sizeof(double));
		Of.write((char *)&v.c[2], sizeof(double));
		Of.close();
	}

	void getSST(Grain *f)
	{
		// Ось четвёртого порядка
		for (double fi = 0; fi < 2 * PI; fi += PI_2)
		{
			saveSSTPoints(f->o, fi, "Polus\\SST001.dat", 0, 0, 1);
			saveSSTPoints(f->o, fi, "Polus\\SST001.dat", 0, 1, 0);
			saveSSTPoints(f->o, fi, "Polus\\SST001.dat", 1, 0, 0);
			saveSSTPoints(f->o, fi, "Polus\\SST001.dat", 0, 0, -1);
			saveSSTPoints(f->o, fi, "Polus\\SST001.dat", 0, -1, 0);
			saveSSTPoints(f->o, fi, "Polus\\SST001.dat", -1, 0, 0);
		}
		// Ось второго порядка
		for (double fi = 0; fi < 2 * PI; fi += PI)
		{
			saveSSTPoints(f->o, fi, "Polus\\SST011.dat", 0, 1, 1);
			saveSSTPoints(f->o, fi, "Polus\\SST011.dat", 0, 1, -1);
			saveSSTPoints(f->o, fi, "Polus\\SST011.dat", 1, 0, 1);
			saveSSTPoints(f->o, fi, "Polus\\SST011.dat", 1, 1, 0);
			saveSSTPoints(f->o, fi, "Polus\\SST011.dat", 0, -1, 1);
			saveSSTPoints(f->o, fi, "Polus\\SST011.dat", -1, 1, 0);
			saveSSTPoints(f->o, fi, "Polus\\SST011.dat", -1, 0, 1);
			saveSSTPoints(f->o, fi, "Polus\\SST011.dat", -1, 0, -1);
			saveSSTPoints(f->o, fi, "Polus\\SST011.dat", 1, 0, -1);
			saveSSTPoints(f->o, fi, "Polus\\SST011.dat", 1, -1, 0);
			saveSSTPoints(f->o, fi, "Polus\\SST011.dat", -1, -1, 0);
			saveSSTPoints(f->o, fi, "Polus\\SST011.dat", 0, -1, -1);
		}
		// Ось третьего порядка
		for (double fi = 0; fi < 2 * PI; fi += 2.0 * PI / 3.0)
		{
			saveSSTPoints(f->o, fi, "Polus\\SST111.dat", 1, 1, 1);
			saveSSTPoints(f->o, fi, "Polus\\SST111.dat", 1, 1, -1);
			saveSSTPoints(f->o, fi, "Polus\\SST111.dat", 1, -1, 1);
			saveSSTPoints(f->o, fi, "Polus\\SST111.dat", -1, 1, 1);
			saveSSTPoints(f->o, fi, "Polus\\SST111.dat", 1, -1, -1);
			saveSSTPoints(f->o, fi, "Polus\\SST111.dat", -1, 1, -1);
			saveSSTPoints(f->o, fi, "Polus\\SST111.dat", -1, -1, 1);
			saveSSTPoints(f->o, fi, "Polus\\SST111.dat", -1, -1, -1);
		}

	}
}