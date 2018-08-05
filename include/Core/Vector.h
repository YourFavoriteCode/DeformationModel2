#pragma once
#ifndef CORE_VECTOR_H
#define CORE_VECTOR_H

#include "CommonCore.h"

/***********************************************************************
****************                  ������               *****************
***********************************************************************/

namespace model
{
	class Vector
	{
	public:
		double C[DIM];					//���������� �������

		double getNorm();				//���������� ����� �������
		int Normalize();				//����������� ������
		void setZero();					//�������� ���������� �������
		void set(double,
			double, double);			//����� �������� ���������
		double ScalMult(Vector v);		//��������� ������������ ��������
		Vector VectMult(Vector v);		//��������� ������������ ��������

		Vector operator + (Vector);		//�������� ��������
		Vector operator - (Vector);		//�������� ���������
		Vector operator - ();			//������� �����
		Vector operator * (double);
		void operator += (Vector);		//�������� ����������� �������
		void operator -= (Vector);		//�������� ��������� �������
		void operator *= (double);		//�������� ��������� ������� �� �����
		int operator /= (double);		//�������� ������� ������� �� �����

		Vector();
		~Vector();

	private:

	};

}

#endif