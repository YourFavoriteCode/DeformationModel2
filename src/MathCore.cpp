// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"
#include <cmath>

#include "CommonCore.h"
#include "Vector.h"
#include "Tensor.h"
#include "Tensor4.h"
#include "SlipSystem.h"
#include "Functions.h"
#include "Eigen/Eigen"

namespace model
{
	/***********************************************************************
	*****************      Псевдо-тензор Леви-Чивита     ******************
	***********************************************************************/
	int LeviCivit(const int i, const int j, const int k)
	{
		if (i == 0 && j == 1 && k == 2 || i == 2 && j == 0 && k == 1 || i == 1 && j == 2 && k == 0) return 1;
		else if (i == 2 && j == 1 && k == 0 || i == 0 && j == 2 && k == 1 || i == 1 && j == 0 && k == 2) return -1;
		else return 0;
	}

	/***********************************************************************
	****************      Методы для работы с вектором     *****************
	***********************************************************************/
	Vector::Vector()
	{
		c[0] = 0;
		c[1] = 0;
		c[2] = 0;
	}

	void Vector::set(const double x, const double y, const double z)
	{
		c[0] = x;
		c[1] = y;
		c[2] = z;
	}

	double Vector::getNorm()
	{
		double sum = 0;
		for (int i = 0; i < DIM; i++)
		{
			sum += c[i] * c[i];
		}
		double res;	
		if (isNormalDouble(sum)) res = sqrt(sum);
		else res = -1;//В случае ошибки вернёт -1
		return res;
	}

	int Vector::normalize()
	{
		double norm = getNorm();
		if (norm > 0)
		{
			for (int i = 0; i < DIM; i++)
			{
				c[i] /= norm;
			}
			return 0;
		}
		else
		{
			return -1;		//Код ошибки деления на 0 и QNaN
		}
		
	}

	void Vector::setZero()
	{
		for (int i = 0; i < DIM; i++)
		{
			c[i] = 0;
		}
	}

	double Vector::scalMult(const Vector v)
	{
		double res = 0;
		for (int i = 0; i < DIM; i++)
		{
			res += c[i] * v.c[i];
		}
		return res;
	}
	
	Vector Vector::vectMult(const Vector v)
	{
		Vector res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					res.c[i] += LeviCivit(i, j, k) * c[j] * v.c[k];
				}
			}
		}
		return res;
	}

	Vector Vector::operator+(const Vector v)
	{
		Vector buf;
		for (int i = 0; i < DIM; i++)
		{
			buf.c[i] = c[i] + v.c[i];
		}
		return buf;
	}

	Vector Vector::operator-(const Vector v)
	{
		Vector buf;
		for (int i = 0; i < DIM; i++)
		{
			buf.c[i] = c[i] - v.c[i];
		}
		return buf;
	}
	
	Vector Vector::operator-()
	{
		Vector buf;
		for (int i = 0; i < DIM; i++)
		{
			buf.c[i] = -c[i];
		}
		return buf;
	}

	Vector Vector::operator*(const double r)
	{
		Vector buf;
		for (int i = 0; i < DIM; i++)
		{
			buf.c[i] = c[i] * r;
		}
		return buf;
	}

	void Vector::operator += (const Vector v)
	{
		for (int i = 0; i < DIM; i++)
		{
			c[i] += v.c[i];
		}
	}
	void Vector::operator -= (const Vector v)
	{
		for (int i = 0; i < DIM; i++)
		{
			c[i] -= v.c[i];
		}
	}
	void Vector::operator *= (const double r)
	{
		for (int i = 0; i < DIM; i++)
		{
			c[i] *= r;
		}
	}

	int Vector::operator /= (const double r)
	{
		if (r != 0)
		{
			for (int i = 0; i < DIM; i++)
			{
				c[i] /= r;
			}
		}
		else
		{
			return -1;
		}
		return 0;
	}

	Vector::~Vector()
	{
	}

	/***********************************************************************
	****************      Методы для работы с тензором     *****************
	***********************************************************************/

	void Tensor::transp()
	{
		double buf;
		for (int i = DIM-1; i >= 0; --i)
		{
			for (int j = i - 1; j >= 0; --j)
			{
				buf = c[i][j];
				c[i][j] = c[j][i];
				c[j][i] = buf;
			}
		}
	}

	double Tensor::getDet()
	{
		return (c[0][0] * c[1][1] * c[2][2] + c[0][1] * c[1][2] * c[2][0] +
			c[1][0] * c[2][1] * c[0][2] - c[2][0] * c[1][1] * c[0][2] -
			c[0][0] * c[1][2] * c[2][1] - c[2][2] * c[0][1] * c[1][0]);
	}

	void Tensor::set(double C00, double C01, double C02,
		double C10, double C11, double C12,
		double C20, double C21, double C22)
	{
		c[0][0] = C00;
		c[0][1] = C01;
		c[0][2] = C02;

		c[1][0] = C10;
		c[1][1] = C11;
		c[1][2] = C12;

		c[2][0] = C20;
		c[2][1] = C21;
		c[2][2] = C22;
	}

	void Tensor::setZero()
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				c[i][j] = 0;
			}
		}
	}

	Vector Tensor::getRow(const int n)
	{
		Vector res;
		res.set(c[n][0], c[n][1], c[n][2]);
		return res;
	}

	Vector Tensor::getCol(const int n)
	{
		Vector res;
		res.set(c[0][n], c[1][n], c[2][n]);
		return res;
	}

	void Tensor::setUnit()
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				c[i][j] = (i == j) ? 1.0 : 0;
			}
		}
	}

	double Tensor::doubleScalMult(Tensor t)
	{
		double res = 0;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				res += c[i][j] * t.c[j][i];
			}
		}
		return res;
	}

	Tensor Tensor::getSymmetryPart()
	{
		Tensor res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				res.c[i][j] = c[i][j] + c[j][i];
			}
		}
		res /= 2.0;
		return res;
	}

	Tensor Tensor::getAntiSymmetryPart()
	{
		Tensor res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				res.c[i][j] = c[i][j] - c[j][i];
			}
		}
		res /= 2.0;
		return res;
	}
	
	void Tensor::rotationMatrix(double fi, Vector v)
	{
		double cosFi = cos(fi);
		double sinFi = sin(fi);
		double cosFi1 = (1 - cosFi);

		c[0][0] = cosFi1 * v.c[0] * v.c[0] + cosFi;
		c[0][1] = cosFi1 * v.c[0] * v.c[1] - sinFi * v.c[2];
		c[0][2] = cosFi1 * v.c[0] * v.c[2] + sinFi * v.c[1];

		c[1][0] = cosFi1 * v.c[1] * v.c[0] + sinFi * v.c[2];
		c[1][1] = cosFi1 * v.c[1] * v.c[1] + cosFi;
		c[1][2] = cosFi1 * v.c[1] * v.c[2] - sinFi * v.c[0];

		c[2][0] = cosFi1 * v.c[2] * v.c[0] - sinFi * v.c[1];
		c[2][1] = cosFi1 * v.c[2] * v.c[1] + sinFi * v.c[0];
		c[2][2] = cosFi1 * v.c[2] * v.c[2] + cosFi;
	}

	void Tensor::getAxisAngle(double* fi, Vector* v)
	{
		double tr = c[0][0] + c[1][1] + c[2][2];
		double theta = acos((tr - 1) * 0.5);

		double omPreCalc = 1.0 / (2 * sin(theta));

		Vector w;
		w.c[0] = omPreCalc * (c[2][1] - c[1][2]);
		w.c[1] = omPreCalc * (c[0][2] - c[2][0]);
		w.c[2] = omPreCalc * (c[1][0] - c[0][1]);

		*fi = theta;
		*v = w;
	}

	Tensor::Tensor()
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				c[i][j] = 0;
			}
		}
	}

	Tensor Tensor::operator + (const Tensor t)
	{
		Tensor res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				res.c[i][j] = c[i][j] + t.c[i][j];
			}
		}
		return res;
	}

	void Tensor::operator += (const Tensor t)
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				c[i][j] += t.c[i][j];
			}
		}
	}

	Tensor Tensor::operator - (const Tensor t)
	{
		Tensor res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				res.c[i][j] = c[i][j] - t.c[i][j];
			}
		}
		return res;
	}

	Tensor Tensor::operator - ()
	{
		Tensor res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				res.c[i][j] = -c[i][j];
			}
		}
		return res;
	}

	void Tensor::operator -= (const Tensor t)
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				c[i][j] -= t.c[i][j];
			}
		}
	}

	Tensor operator *(const double r, Tensor &t)
	{
		Tensor res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				res.c[i][j] = t.c[i][j]*r;
			}
		}
		return res;
	}

	Tensor operator *(Tensor &t, const double r)
	{
		Tensor res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				res.c[i][j] = t.c[i][j] * r;
			}
		}
		return res;
	}

	Tensor Tensor::operator *(const Tensor t)
	{
		Tensor res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					res.c[i][j] += c[i][k] * t.c[k][j];
				}
			}
		}
		return res;
	}

	void Tensor::operator *=(const Tensor t)
	{
		Tensor buf;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					buf.c[i][j] += c[i][k] * t.c[k][j];
				}
				c[i][j] = buf.c[i][j];
			}
		}
	}

	void Tensor::operator *= (const double r)
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				c[i][j] *= r;
			}
		}
	}

	int Tensor::operator /= (const double r)
	{
		if (r != 0)
		{
			for (int i = 0; i < DIM; i++)
			{
				for (int j = 0; j < DIM; j++)
				{
					c[i][j] /= r;
				}
			}
		}
		else
		{
			return -1;	//Код ошибки деления на 0
		}
		return 0;
	}

	Tensor::~Tensor()
	{
	}

	/***********************************************************************
	***********      Методы для работы с системой скольжения     ***********
	***********************************************************************/

	void SlipSystem::Initialize(int nx, int ny,
		int nz, int bx, int by, int bz)
	{
		/*
		Инициализация трёхиндексовых СС
		*/
		n.set(nx, ny, nz);
		b.set(bx, by, bz);
		
		n.normalize();
		b.normalize();
	}

	void SlipSystem::Initialize(int nx1, int nx2, int nx3, int nx4,
		int bx1, int bx2, int bx3, int bx4, int c)
	{
		/*
		Инициализация четырёхиндексовых СС
		*/
		n.set(nx1, nx2, nx4);

		double u = (2.0 * bx1 + bx2) / 3.0;
		double v = (bx1 + 2.0 * bx2) / 3.0;
		double w = bx4 / 3.0;
		b.set(u, v, w);

		n.normalize();
		b.normalize();
	}

	SlipSystem::SlipSystem()
	{
		t = 0;
		tc = 0;
		tbs = 0;
		dgm = 0;
		gmm = 0;
	}

	SlipSystem::~SlipSystem()
	{

	}

	/***********************************************************************
	********      Методы для работы с тензором четвёртого ранга    *********
	***********************************************************************/
	
	Tensor4::Tensor4()
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					for (int l = 0; l < DIM; l++)
					{
						C[i][j][k][l] = 0;
					}
				}
			}
		}

	}

	Tensor4::~Tensor4()
	{

	}
	

	void Tensor4::setZero()
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					for (int l = 0; l < DIM; l++)
					{
						C[i][j][k][l] = 0;
					}
				}
			}
		}
	}

	void Tensor4::Symmetrize()
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					for (int l = 0; l < DIM; l++)
					{
						double buf = (C[i][j][k][l] + C[j][i][k][l] + C[i][j][l][k] + C[j][i][l][k]) / 4.0;
						C[i][j][k][l] = C[j][i][k][l] = C[i][j][l][k] = C[j][i][l][k] = buf;

					}
				}
			}
		}
	}
	
	Tensor4 Tensor4::ToLSK(Tensor O)
	{
		Tensor4 e;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					for (int l = 0; l < DIM; l++)
					{
						for (int i1 = 0; i1 < DIM; i1++)
						{
							for (int j1 = 0; j1 < DIM; j1++)
							{
								for (int k1 = 0; k1 < DIM; k1++)
								{
									for (int l1 = 0; l1 < DIM; l1++)
									{
										e.C[i][j][k][l] += O.c[i][i1] * O.c[j][j1] * O.c[k][k1] * O.c[l][l1] * C[i1][j1][k1][l1];
									}
								}
							}
						}
					}
				}
			}
		}
		return e;
	}

	void Tensor4::operator += (Tensor4 e)
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					for (int l = 0; l < DIM; l++)
					{
						C[i][j][k][l] += e.C[i][j][k][l];
					}
				}
			}
		}
	}

	void Tensor4::operator -= (Tensor4 e)
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					for (int l = 0; l < DIM; l++)
					{
						C[i][j][k][l] -= e.C[i][j][k][l];
					}
				}
			}
		}
	}

	void Tensor4::operator *= (double r)
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					for (int l = 0; l < DIM; l++)
					{
						C[i][j][k][l] *= r;
					}
				}
			}
		}
	}

	int Tensor4::operator /= (double r)
	{
		if (r != 0)
		{
			for (int i = 0; i < DIM; i++)
			{
				for (int j = 0; j < DIM; j++)
				{
					for (int k = 0; k < DIM; k++)
					{
						for (int l = 0; l < DIM; l++)
						{
							C[i][j][k][l] /= r;
						}
					}
				}
			}
		}
		else
		{
			return -1;
		}
		return 0;
	}

	/***********************************************************************
	***********   Прочие функции, не вошедшие в состав классов   ***********
	***********************************************************************/
	
	Tensor VectMult(Vector v, Tensor t)
	{
		Tensor res;
		for (int k = 0; k < DIM; k++)
		{
			Vector buf = t.getCol(k);
			Vector r = v.vectMult(buf);
			for (int i = 0; i < DIM; i++)
			{
				res.c[i][k] = r.c[i];
			}
		}
		return res;
	}
	
	Tensor VectMult(Tensor t, const Vector v)
	{
		Tensor res;
		for (int k = 0; k < DIM; k++)
		{
			Vector buf = t.getRow(k);
			Vector r = buf.vectMult(v);
			for (int i = 0; i < DIM; i++)
			{
				res.c[k][i] = r.c[i];
			}
		}
		return res;
	}

	Vector ScalMult(Vector v, Tensor t)
	{
		Vector res;
		for (int k = 0; k < DIM; k++)
		{
			Vector buf = t.getCol(k);
			res.c[k] = v.scalMult(buf);
		}
		return res;
	}

	Vector ScalMult(Tensor t, const Vector v)
	{
		Vector res;
		for (int k = 0; k < DIM; k++)
		{
			Vector buf = t.getRow(k);
			res.c[k] = buf.scalMult(v);
		}
		return res;
	}

	Tensor Transp(Tensor t)
	{
		Tensor res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				res.c[i][j] = t.c[j][i];
			}
		}
		return res;
	}

	double Tensor::getL(int n)
	{
		Eigen::MatrixXf m(DIM, DIM);
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				m(i, j) = c[i][j];
			}
		}
		Eigen::EigenSolver<Eigen::MatrixXf> solver(m);
		return solver.eigenvalues().col(0)[n].real();
	}
}