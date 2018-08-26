#pragma once
#ifndef __GRAIN_STRUCTURE_H 
#define __GRAIN_STRUCTURE_H

#include "Polycrystall.h"

namespace model
{

	class GrainStructure
	{
	public:
		
		std::vector<Vector> centerPos;

		// Создает простую трехмерную структуру из зерен-кубиков 
		void makeCubicStructure(Polycrystall*);
		// Создает трехмерную структуру из многогранников Вороного
		void makeVoronoiStructure(Polycrystall*);
	private:
	};
}
#endif