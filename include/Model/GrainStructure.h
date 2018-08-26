#pragma once
#ifndef __GRAIN_STRUCTURE_H 
#define __GRAIN_STRUCTURE_H

#include "Polycrystall.h"

namespace model
{

	class GrainStructure
	{
	public:
		void makeCubicStruct(Polycrystall*);
		void makeVoronoiStructure(Polycrystall*);
	private:
	};
}
#endif