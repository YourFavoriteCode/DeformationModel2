// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"

#include "SlipSystem.h"

namespace model
{
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
}