// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"

#include "Tensor4.h"

namespace model
{
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
						c[i][j][k][l] = 0;
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
						c[i][j][k][l] = 0;
					}
				}
			}
		}
	}

	void Tensor4::symmetrize()
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					for (int l = 0; l < DIM; l++)
					{
						double buf = (c[i][j][k][l] + c[j][i][k][l] + c[i][j][l][k] + c[j][i][l][k]) / 4.0;
						c[i][j][k][l] = c[j][i][k][l] = c[i][j][l][k] = c[j][i][l][k] = buf;

					}
				}
			}
		}
	}

	Tensor4 Tensor4::toLsk(Tensor O)
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
										e.c[i][j][k][l] += O.c[i][i1] * O.c[j][j1] * O.c[k][k1] * O.c[l][l1] * c[i1][j1][k1][l1];
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
						c[i][j][k][l] += e.c[i][j][k][l];
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
						c[i][j][k][l] -= e.c[i][j][k][l];
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
						c[i][j][k][l] *= r;
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
							c[i][j][k][l] /= r;
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

	Tensor Tensor4::doubleScalMult(Tensor t)
	{
		Tensor res;
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < DIM; k++)
				{
					for (int l = 0; l < DIM; l++)
					{
					res.c[i][j] += c[i][j][k][l] * t.c[l][k];
					}
				}
			}
		}
		return res;
	}
}