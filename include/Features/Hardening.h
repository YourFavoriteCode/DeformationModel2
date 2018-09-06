#pragma once
#ifndef __HARDENING_H 
#define __HARDENING_H

#include "Grain.h"
/*
*Механизмы упрочнения
*/
namespace model
{

	void Base_hardening(Grain*);			//Базовое слагаемое упрочнения
	void Boundary_hardening(Grain*);		//Зернограничное упрочнение
	
}
#endif __HARDENING_H
