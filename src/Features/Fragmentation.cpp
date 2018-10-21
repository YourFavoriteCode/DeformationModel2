// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"
#include "GrainStructure.h"
#include "Distributions.h"

namespace model
{
	void calcDislocationsPower()
	{

	}
	

	// Критерий Треска-Геста (максимального касательного напряжения)
	bool criteria(Grain g)
	{
		double maxStress = prms::fragmentationCriteria;
		double s[3];
		for (int i = 0; i < 3; i++)
		{
			s[i] = g.sgm.getL(i);
		}
		double t0 = (s[1] - s[0]) / 2;
		double t1 = (s[2] - s[0]) / 2;
		double t2 = (s[2] - s[1]) / 2;
		return (abs(t0) >= maxStress || abs(t1) >= maxStress || abs(t2) >= maxStress);
	}

	double getRnd() { return double(rand()) / RAND_MAX; }

	void GrainStructure::fragmentate()
	{
		// Модификация массива в цикле - это плохо
		for (Grain g : polycrystall->c)
		{
			if (criteria(g))
			{
				split(g.position);
				return;
			}
		}
	}

	void deleteIndex(std::map<int, Vector> *map, int num)
	{
		std::map<int, Vector>::iterator it;
		it = map->find(num);
		map->erase(it);	
	}

	void GrainStructure::split(int oldPos)
	{
		//TODO: тут ищется плоскость разделения

		// Удаляем номера и центры из индекса
		deleteIndex(&posMap, oldPos);

		// Этот рандом чисто дял тестов
		double x = xMin + getRnd()*(xMax - xMin);
		double y = yMin + getRnd()*(yMax - yMin);
		double z = zMin + getRnd()*(zMax - zMin);
		Vector vector1(x, y, z);
		x = xMin + getRnd()*(xMax - xMin);
		y = yMin + getRnd()*(yMax - yMin);
		z = zMin + getRnd()*(zMax - zMin);
		Vector vector2(x, y, z);

		// Для оптимизации старое зерно вместо удаления будет вторым новым субзерном
		int n1 = lastId++, n2 = lastId++;
		int npos = polycrystall->findGrainPosById(oldPos);
		polycrystall->c[npos].position = n1;
		// Новое субзерно полностью сохранит состояние старого
		// Изменится лишь его идентификатор, ориентация решетки и форма 
		Grain newGrain = polycrystall->c[npos];
		newGrain.position = n2;
		polycrystall->c.push_back(newGrain);
		polycrystall->totalGrainCount++;
		posMap.insert(std::pair<int, Vector>(n1, Vector()));
		posMap.insert(std::pair<int, Vector>(n2, Vector()));
		updateContainer();

	}
}