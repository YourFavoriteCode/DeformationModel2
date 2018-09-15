#pragma once
#ifndef __GRAIN_STRUCTURE_H 
#define __GRAIN_STRUCTURE_H

#include "Polycrystall.h"
#include "voro++.hh"

namespace model
{

	class GrainStructure
	{
	public:
		
		GrainStructure(Polycrystall*);

		std::vector<Vector> centerPos;

		// Создает трехмерную структуру из многогранников Вороного
		void makeVoronoiStructure();

		void updateStructure(voro::container*);
	private:
		Polycrystall *polycrystall;
	};
}
#endif