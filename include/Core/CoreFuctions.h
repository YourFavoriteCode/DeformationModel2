#pragma once
#ifndef CORE_FUNCTIONS_H
#define CORE_FUNCTIONS_H

#include "Vector.h"
#include "Tensor.h"

namespace model
{
	int LeviCivit(int i, int j, int k);	//������-������ ����-������

	Tensor VectMult(Vector, Tensor);	//��������� ������������ ������� �� ������
	Tensor VectMult(Tensor, Vector);	//��������� ������������ ������� �� ������
	Vector ScalMult(Vector, Tensor);	//��������� ������������ ������� �� ������
	Vector ScalMult(Tensor, Vector);	//��������� ������������ ������� �� ������

	Tensor Transp(Tensor);				//���������� ����������������� ��������
}

#endif