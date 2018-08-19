// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"

#include "cmath"

#include "Vector.h"
#include "Functions.h"

namespace model
{
	/***********************************************************************
	****************      Методы для работы с вектором     *****************
	***********************************************************************/

	Vector::Vector()
	{
		c[0] = 0;
		c[1] = 0;
		c[2] = 0;
	}

	Vector::Vector(const double x, const double y, const double z)
	{
		c[0] = x;
		c[1] = y;
		c[2] = z;
	}

	Vector::~Vector()
	{
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
		return isNormalDouble(sum) ? sqrt(sum) : -1;
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

}