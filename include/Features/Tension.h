#pragma once
#ifndef __TENSION_H 
#define __TENSION_H

#include "Tensor.h"
#include "Tensor4.h"

/*
*Функции для опредления тензоров напряжений и деформаций при одноосном растяжении
*/

namespace model
{
	Tensor TensionStrainCalc(Tensor4 P, 
		Tensor D_in, double tens_comp);					//Вычисление тензора деформаций при разгрузке
	Tensor TensionStressCalc(Tensor4 &P,
		Tensor &D_in, Tensor &D);						//Вычисление тензора напряжений при растяжении
	Tensor UnloadingStrainCalc(Tensor4 &P, 
		Tensor &D_in, Tensor &Sgm, double lam);			//Вычисляет тензор деформации при разгрузке
}
#endif __TENSION_H