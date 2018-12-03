// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"
#include "GrainStructure.h"
#include "Distributions.h"
#include "Rotations.h"

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
		double t0 = abs(s[1] - s[0]) / 2;
		double t1 = abs(s[2] - s[0]) / 2;
		double t2 = abs(s[2] - s[1]) / 2;
		return (t0 >= maxStress || t1 >= maxStress || t2 >= maxStress);
	}

	double getRnd() { return double(rand()) / RAND_MAX; }

	void GrainStructure::fragmentate()
	{
		// Случайным образом выбирается стартовый элемент и направление обхода 
		std::vector<Grain>::iterator iter;
		int randStart = std::rand() % polycrystall->c.size();
		bool reverse = std::rand() > 0.5;
		for (iter = polycrystall->c.begin()+randStart;
			iter != (reverse ? polycrystall->c.begin() : polycrystall->c.end());
			reverse ? iter-- : iter++)
		{
			if (criteria(*iter))
			{
				split((*iter).position);
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

	Vector getIndexVector(std::map<int, Vector> *map, int num)
	{
		std::map<int, Vector>::iterator it;
		it = map->find(num);
		return it->second;
	}

	/*
	Методом перебора ищет плоскость максимального касательного
	напряжения, возвращает нормаль к этой плоскости
	*/
	Vector maxStressPlane(Tensor sgm)
	{
		Vector nMax;
		double mMax = 0;

		//Генерируем случайные плоскости
		for (int i = 0; i < 200; i++)
		{
			//Случайное направление
			double a = -1 + ((double)rand() / RAND_MAX) * 2;
			double b = -1 + ((double)rand() / RAND_MAX) * 2;
			double c = -1 + ((double)rand() / RAND_MAX) * 2;
			Vector tau(a,b,c);				//Касательный вектор
			double ksi = -(-a + b) / c;		//Гарантируем ортогональность
			Vector n(1.0, 1.0, ksi);		//Вектор нормали
			tau.normalize();
			n.normalize();
			double m = scalMult(tau, sgm).scalMult(n);
			if (abs(m) > mMax)
			{
				nMax = n;
			}
		}
		return nMax;
	}

	/*
	Определяет позиции новых центров многогранников
	и добавляет их в общую мезо-структуру
	*/
	void GrainStructure::split(int oldPos)
	{
		// Для оптимизации старое зерно вместо удаления будет вторым новым субзерном
		int n1 = lastId++, n2 = lastId++;
		int npos = polycrystall->findGrainPosById(oldPos);
		polycrystall->c[npos].position = n1;
		polycrystall->c[npos].sgm /= 10;
		rotateByTaylor(&(polycrystall->c[npos]));
		// Поиск плоскости разделения
		Vector nSplit = maxStressPlane(polycrystall->c[npos].sgm);
		// Предполагаем, что плоскость проходит через
		// центр многогранника Вороного
		double coef = 1e-4;
		Vector oldVec = getIndexVector(&posMap, oldPos);
		// Новые центры при этом будут смещены от старого
		// в обе стороны по линии нормали на малую величину 
		Vector newVec1 = oldVec - nSplit * coef;
		Vector newVec2 = oldVec + nSplit * coef;
		// Удаляем номер и центр из индекса
		deleteIndex(&posMap, oldPos);
			
		// Новое субзерно полностью сохранит состояние старого
		// Изменится лишь его идентификатор, ориентация решетки и форма 
		Grain newGrain = polycrystall->c[npos];
		newGrain.position = n2;
		// Нужно поворачивать и новое и старое зерно
		rotateByTaylor(&newGrain);
		polycrystall->c.push_back(newGrain);
		polycrystall->totalGrainCount++;
		posMap.insert(std::pair<int, Vector>(n1, newVec1));
		posMap.insert(std::pair<int, Vector>(n2, newVec2));
		updateContainer();
	}
}