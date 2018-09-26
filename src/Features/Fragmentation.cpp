// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"
#include "GrainStructure.h"

namespace model
{
	void calcDislocationsPower()
	{

	}
	
	void GrainStructure::fragmentate()
	{
		for (Grain g : polycrystall->c)
		{
			if (true)
			{
				split(&g);
			}
		}
	}

	double getRnd() { return double(rand()) / RAND_MAX; }

	void deleteIndex(std::map<int, Vector> map, int num)
	{
		std::map<int, Vector>::iterator it;
		it = map.find(num);
		map.erase(it);	
	}

	void GrainStructure::split(Grain* oldGrain)
	{
		//TODO: тут ищется плоскость разделения

		// Удаляем номера и центры из индекса
		int id = oldGrain->position;
		deleteIndex(posMap, id);

		// Этот рандом чисто дял тестов
		double x = xMin + getRnd()*(xMax - xMin);
		double y = yMin + getRnd()*(yMax - yMin);
		double z = zMin + getRnd()*(zMax - zMin);
		Vector vector1(x, y, z);
		x = xMin + getRnd()*(xMax - xMin);
		y = yMin + getRnd()*(yMax - yMin);
		z = zMin + getRnd()*(zMax - zMin);
		Vector vector2(x, y, z);

		// Новое субзерно полностью сохранит состояние старого
		// Изменится лишь его идентификатор, ориентация решетки и форма 
		Grain newGrain = *oldGrain; 
		// Для оптимизации старое зерно вместо удаления будет вторым новым субзерном
		int n1 = lastId++, n2 = lastId++;
		newGrain.position = n1;
		(*oldGrain).position = n2;
		polycrystall->c.push_back(newGrain);
		polycrystall->totalGrainCount++;
		posMap.insert(std::pair<int, Vector>(n1, Vector()));
		posMap.insert(std::pair<int, Vector>(n2, Vector()));
		updateContainer();


	}
}