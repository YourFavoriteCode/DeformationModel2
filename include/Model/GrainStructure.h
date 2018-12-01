#pragma once
#ifndef __GRAIN_STRUCTURE_H 
#define __GRAIN_STRUCTURE_H

#include "Polycrystall.h"
#include "voro++.hh"
#include "map"

namespace model
{

	class GrainStructure
	{
	public:
		
		GrainStructure(Polycrystall*, bool);

		voro::container* container;
		std::map<int, Vector> posMap;
		bool periodic;

		// Создает трехмерную структуру из многогранников Вороного
		void makeVoronoiStructure();

		void updateStructure(voro::container*);

		void printStructureInfo(voro::container *, const char*);

		void updateContainer();

		void fragmentate();

	private:
		Polycrystall *polycrystall;		// Ссылка на поликристалл
		int lastId;						// Последний зарегистрированный идентификатор
		void split(int);
		voro::container* makeContainer();
		// Геометрия контейнера для расчетной области
		const double xMin = 0, xMax = 1;
		const double yMin = 0, yMax = 1;
		const double zMin = 0, zMax = 1;
	};
}
#endif