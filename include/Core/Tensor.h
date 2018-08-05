#pragma once
#ifndef CORE_TENSOR_H
#define CORE_TENSOR_H

#include "CommonCore.h"
#include "Vector.h"

/***********************************************************************
****************                  ������               *****************
***********************************************************************/

namespace model
{

	class Tensor
	{
	public:
		double C[DIM][DIM];				//���������� �������

		double getDet();				//���������� ������������ ������� ���������
		void set(double, double, double,
			double, double, double,
			double, double, double);	//������� ���� ��������� �������
		void setZero();					//�������� ���������� �������
		void setUnit();					//������ ������� ��������� ������� ���������
		void Transp();					//������������� ������� ��������� �������
		double doubleScalMult(Tensor);	//������ (������� ��������� ������������ ��������)
		void RotMatr(double, Vector);	//������������� ������ �������� ������ �������� ��� �� �������� ����
		void getAxisAngle(double*, Vector*);	//���������� ��� � ���� ��������, �� ������� ��� ��������� ������

		Tensor getSymmetryPart();		//���������� ������������ ����� �������
		Tensor getAntiSymmetryPart();	//���������� ���������������� ����� �������
		Vector getRow(int);				//������ �� ��������� �������� ������
		Vector getCol(int);				//������ �� ��������� ��������� ������

		double getL(int n);				//���������� ����������� ����� � ������� n

		Tensor operator + (Tensor);		//�������� �������� ��������
		Tensor operator - (Tensor);		//�������� ��������� ��������
		Tensor operator - ();			//������� �����
		Tensor operator * (Tensor);		//�������� ��������� ��������

		friend Tensor operator * (Tensor&, double);		//�������� ��������� �� �����
		friend Tensor operator * (double, Tensor&);		//�������������

		void operator += (Tensor);		//�������� ����������� �������
		void operator -= (Tensor);		//�������� ��������� �������
		void operator *= (Tensor);		//�������� ���������� �� ������
		void operator *= (double);		//�������� ��������� ������� �� �����
		int operator /= (double);		//�������� ������� ������� �� �����

		Tensor();
		~Tensor();

	private:

	};
}

#endif