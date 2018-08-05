#pragma once
#ifndef __HARDENING_H 
#define __HARDENING_H

#include "Fragment.h"
/*
*Механизмы упрочнения
*/
namespace model
{

	void Base_hardening(Fragment*);			//Базовое слагаемое упрочнения
	void Boundary_hardening(Fragment*);		//Зернограничное упрочнение
	
}
#endif __HARDENING_H
