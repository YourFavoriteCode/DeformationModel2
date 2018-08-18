#include "stdafx.h"

#include "Loading.h"

namespace model {

	Loading::Loading()
	{
		cycle = 0;
		cycleCount = 1;
		currentStep = 0;
		tensionComponent = 0;
		maxStress = 1e4;
		paramLambda = 2.5;
		additionalStrain = 1e-5;
	}

	Loading::~Loading()
	{

	}

	void Loading::setLoad(Tensor t, int count) {
		gradV = t;
		cycleCount = count;
	}
}