// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"

#include "cmath"

#include "Tensor.h"
#include "Eigen/Eigen"

namespace model
{
	/***********************************************************************
	****************      Методы для работы с тензором     *****************
	***********************************************************************/

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

	Tensor::~Tensor()
	{
	}

	void Tensor::transp()
	{
		double buf;
		for (int i = DIM - 1; i >= 0; --i)
		{
			for (int j = i - 1; j >= 0; --j)
			{
				buf = c[i][j];
				c[i][j] = c[j][i];
				c[j][i] = buf;
			}
		}
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

	double Tensor::getDet()
	{
		return (c[0][0] * c[1][1] * c[2][2] + c[0][1] * c[1][2] * c[2][0] +
			c[1][0] * c[2][1] * c[0][2] - c[2][0] * c[1][1] * c[0][2] -
			c[0][0] * c[1][2] * c[2][1] - c[2][2] * c[0][1] * c[1][0]);
	}

	void Tensor::set(double c00, double c01, double c02,
		double c10, double c11, double c12,
		double c20, double c21, double c22)
	{
		c[0][0] = c00;
		c[0][1] = c01;
		c[0][2] = c02;

		c[1][0] = c10;
		c[1][1] = c11;
		c[1][2] = c12;

		c[2][0] = c20;
		c[2][1] = c21;
		c[2][2] = c22;
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
		return Vector(c[n][0], c[n][1], c[n][2]);
	}

	Vector Tensor::getCol(const int n)
	{
		return Vector(c[0][n], c[1][n], c[2][n]);
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
				res.c[i][j] = t.c[i][j] * r;
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

	Vector ScalMult(Tensor t, Vector v)
	{
		Vector res;
		for (int k = 0; k < DIM; k++)
		{
			Vector buf = t.getRow(k);
			res.c[k] = buf.scalMult(v);
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

	Tensor VectMult(Tensor t, Vector v)
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

}