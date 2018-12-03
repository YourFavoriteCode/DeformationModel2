// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"
#include <cmath>

#include "Params.h"
#include "Materials.h"
#include "Hardening.h"

namespace model
{
	void hardeningBase(Grain *f)
	{
		double nsd = 0;	//Накопленный сдвиг
		for (int k = 0; k < f->ssCount; k++)
		{
			nsd += f->ss[k].gmm;		//Суммируем по всем СС
		}
		for (int k = 0; k < f->ssCount; k++)
		{
			double osn = 0;
			for (int j = 0; j < f->ssCount; j++)
			{
				double a = (j == k ? prms::hardeningParamBaseA : 1.25 * prms::hardeningParamBaseA);
				if ((f->ss[j].dgm > EPS) && (fabs(nsd) > EPS))
				{
					osn += a*pow(f->ss[j].gmm, prms::hardeningParamBasePsi)*pow((f->ss[j].dgm / prms::shearRateLawDgm0), prms::hardeningParamBaseDelta) / nsd;
				}
			}
			double napr = 0;
			switch (f->materialType)
			{
			case 0://Сталь 45 (ОЦК)
			{
				if (k < 24) napr = STEEL_TC1;
				else if (k < 48) napr = STEEL_TC2;
				else napr = STEEL_TC3;
				break;
			}
			case 1://Медь (ГЦК)
			{
				napr = CUPR_TC;
				break;
			}
			case 2://Титан (ГПУ)
			{
				if (k < 6) napr = TITAN_TC1;
				else if (k < 12) napr = TITAN_TC2;
				else napr = TITAN_TC3;
				break;
			}

			}
			f->ss[k].tc += napr*osn*prms::dt;
		}
	}
	
	void hardeningBoundary(Grain *f)
	{		
		double G = 45.5e+9, nu = 0.35;
		double K, K1, Rmin=5e-7;
		Tensor S(3.333, 0.0, 0.0, 0.0, -0.667, 0.0, 0.0, 0.0, nu*2.667);
		// Цикл по СС текущего фрагмента
		for (int k = 0; k < f->ssCount; k++)	
		{
			if (fabs(f->ss[k].dgm) < EPS) continue;
			// Перевели вектор b текущей СС данного зерна в ЛСК
			Vector b1 = scalMult(f->o, f->ss[k].b);
			double zgu = 0;
			double tbs = 0;
			for (int h = 0; h <	f->neighbors.size(); h++)	//Цикл по фасеткам			
			{
				if (f->ss[k].b.scalMult(f->normals[h]) < 0) continue; //Скольжение от границы - пропускаем
				//double zguk = prms::HARD_BOUND_K * f->ss[k].dgm * f->ss[k].gmm / f->size;
				double min = 1.0;//Минимум
				//min = f->disorientMeasure(h);
				// Цикл по системам соседнего зерна
				for (int p = 0; p < f->neighbors[h]->ssCount; p++)	
				{
					// Перевели вектор b p-ой СС соседнего зерна в ЛСК
					Vector b2 = scalMult(f->neighbors[h]->o, f->neighbors[h]->ss[p].b);
					Vector diff = b1 - b2;
					diff.normalize();
					double M = fabs(diff.scalMult(f->normals[h]));

					if (M < min) min = M;

				}
				K = 0.25;//min/**G*/ / (2 * PI*(1.0 - nu));
				K1 = 30000;// log(f->size / Rmin) / (PI*(f->size - Rmin));
				Tensor OT = f->o;
				OT.transp();
				Vector s11 = scalMult(OT, f->ss[k].n);
				Vector s22 = scalMult(OT, f->ss[k].b);
				s11.normalize();
				s22.normalize();
				tbs += S.doubleScalMult(diadMult(s11, s22));
			}
			tbs *= K*K1;
			f->ss[k].tbs += tbs * prms::hardeningParamBoundK;
		}
		
	}

}