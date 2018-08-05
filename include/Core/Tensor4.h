#pragma once
#ifndef CORE_TENSOR4_H
#define CORE_TENSOR4_H

#include "CommonCore.h"
#include "Tensor.h"

/***********************************************************************
**************         ������ ��������� �����            **************
***********************************************************************/

namespace model
{

	class Tensor4
	{
	public:
		double C[DIM][DIM][DIM][DIM];	//���������� �������

		void setZero();					//��������� ��������� �������
		void Symmetrize();				//������������� ��������� �������

		Tensor4 ToLSK(Tensor O);		//������� ��������� ������� � ���

		void operator += (Tensor4);		//�������� ����������� �������
		void operator -= (Tensor4);		//�������� ��������� �������
		void operator *= (double);		//�������� ��������� ������� �� �����
		int operator /= (double);		//�������� ������� ������� �� �����

		Tensor4();
		~Tensor4();
	private:
	};
}

#endif