// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"

#include "CommonCore.h"

namespace model
{
	//Псевдо-тензор Леви-Чивита третьего ранга
	int LeviCivit(const int i, const int j, const int k)
	{
		if (i == 0 && j == 1 && k == 2 || i == 2 && j == 0 && k == 1 || i == 1 && j == 2 && k == 0) return 1;
		else if (i == 2 && j == 1 && k == 0 || i == 0 && j == 2 && k == 1 || i == 1 && j == 0 && k == 2) return -1;
		else return 0;
	}
}