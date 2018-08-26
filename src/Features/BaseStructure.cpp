// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"

#include "Params.h"
#include "GrainStructure.h"

namespace model
{

	// Псевдо-трехмерная раскладка зерен в виде кубиков с малыми случайными отклонениями
	// нормалей и с периодическими граничными условиями. В зависимости от параметра prms::grainSurroundCount
	// учитываются 6, 14 или 26 близлежащих зерен с различными площадями фасеток контакта.

	void GrainStructure::makeCubicStructure(Polycrystall* poly)
	{
		for (int q = 0; q < poly->totalGrainCount; q++)
		{
			int q1, q2, q3, y;
			get3DPos(q, &q1, &q2, &q3);
			for (int h = 0; h < prms::grainSurroundCount; h++)
			{
				int qq1 = q1, qq2 = q2, qq3 = q3;
				//qq1, qq2, qq3 - координаты зерна-соседа
				//y - номер нормали в соседнем зерне в направлении данного зерна
				double fi = ((double)rand() / RAND_MAX) * (PI / 12);//Случайный угол отклонения нормали
				//TODO: предвычислить наиболее распространенные слагаемые для удобства чтения
				switch (h)
				{
				case 0://Вверх
				{
					poly->c[q].normals[h].set(-sin(fi), sin(fi) / cos(fi), 1 / cos(fi));
					qq3 = q3 == poly->grainCount - 1 ? 0 : q3 + 1;
					y = 5;
					break;
				}
				case 1://От нас
				{
					poly->c[q].normals[h].set(-1 / cos(fi), sin(fi), sin(fi) / cos(fi));
					qq1 = q1 == 0 ? poly->grainCount - 1 : q1 - 1;
					y = 3;
					break;
				}
				case 2://Вправо
				{
					poly->c[q].normals[h].set(sin(fi) / cos(fi), 1 / cos(fi), -sin(fi));
					qq2 = q2 == poly->grainCount - 1 ? 0 : q2 + 1;
					y = 4;
					break;
				}
				case 3://На нас
				{
					poly->c[q].normals[h].set(1 / cos(fi), -sin(fi), sin(fi) / cos(fi));
					qq1 = q1 == poly->grainCount - 1 ? 0 : q1 + 1;
					y = 1;
					break;
				}
				case 4://Влево
				{
					poly->c[q].normals[h].set(sin(fi), -1 / cos(fi), -sin(fi) / cos(fi));
					qq2 = q2 == 0 ? poly->grainCount - 1 : q2 - 1;
					y = 2;
					break;
				}
				case 5://Вниз
				{
					poly->c[q].normals[h].set(sin(fi) / cos(fi), sin(fi), -1 / cos(fi));
					qq3 = q3 == 0 ? poly->grainCount - 1 : q3 - 1;
					y = 0;
					break;
				}
				}
				int index = get1DPos(qq1, qq2, qq3);
				poly->c[q].neighbors[h] = &poly->c[index];					// Устанавливаем ссылки на соседние зерна
				poly->c[index].neighbors[y] = &poly->c[q];					// И в обратную сторону
				poly->c[q].normals[h].normalize();
				poly->c[index].normals[y] = -poly->c[q].normals[h];			// Сохранение нормалей граничащих фасеток
			}
		}
	}
}