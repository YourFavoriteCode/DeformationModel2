#pragma once
#ifndef __HARDENING_H 
#define __HARDENING_H

#include "Grain.h"
/*
*Механизмы упрочнения
*/
namespace model
{

	void hardeningBase(Grain*);			// Базовое слагаемое упрочнения
	void hardeningBoundary(Grain*);		// Зернограничное слагаемое упрочнения
	
}
#endif __HARDENING_H
