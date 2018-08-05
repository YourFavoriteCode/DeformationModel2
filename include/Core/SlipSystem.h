#pragma once
#ifndef SLIP_SYSTEM_H
#define SLIP_SYSTEM_H

#include "Vector.h"

/***********************************************************************
****************         ������� ����������            *****************
***********************************************************************/

namespace model
{

	class SlipSystem
	{
	public:
		Vector n;						//������ �������
		Vector b;						//������ ��������
		double t;						//����������� ����������� ����������
		double tc;						//����������� ����������� ����������
		double dgm;						//�������� ������
		double gmm;						//����������� �����
		double tbs;						//�������� ����������� ����������
		void Initialize(int, int, int,
			int, int, int);	//������������� ��������
		void Initialize(int, int, int, int,
			int, int, int, int, int);

		SlipSystem();
		~SlipSystem();
	};

}

#endif