// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"

#include "Polycrystall.h"
#include "Params.h"

namespace model
{

	// Псевдо-трехмерная раскладка зерен в виде кубиков с малыми случайными отклонениями
	// нормалей и с периодическими граничными условиями. В зависимости от параметра prms::grainSurroundCount
	// учитываются 6, 14 или 26 близлежащих зерен с различными площадями фасеток контакта.

	void Polycrystall::makeGrainStruct()
	{
		for (int q = 0; q < totalGrainCount; q++)
		{
			int q1, q2, q3, y;
			get3DPos(q, &q1, &q2, &q3);
			for (int h = 0; h < prms::grainSurroundCount; h++)
			{
				//Если контакт уже был задан - пропускаем
				if (c[q].contact[h] != -1) continue;
				//Определяем, граничат ли фрагменты
				//Первые 6, т.е. боковые грани, граничат всегда
				double a = h < 6 ? 1 : ((double)rand() / RAND_MAX);//На всё воля божья
				if (a < 0.5)
				{
					//Контакта нет - тоже пропускаем
					c[q].contact[h] = 0;
					continue;
				}
				int qq1 = q1, qq2 = q2, qq3 = q3;
				//qq1, qq2, qq3 - координаты зерна соседа
				//y - номер нормали в соседнем зерне в направлении данного зерна
				double fi = ((double)rand() / RAND_MAX) * (PI / 12);//Случайный угол отклонения нормали
																	//TODO: предвычислить наиболее распространенные слагаемые для удобства чтения
				switch (h)
				{
				case 0://Вверх
				{
					c[q].normals[h].set(-sin(fi), sin(fi) / cos(fi), 1 / cos(fi));
					qq3 = q3 == grainCount - 1 ? 0 : q3 + 1;
					y = 5;
					break;
				}
				case 1://От нас
				{
					c[q].normals[h].set(-1 / cos(fi), sin(fi), sin(fi) / cos(fi));
					qq1 = q1 == 0 ? grainCount - 1 : q1 - 1;
					y = 3;
					break;
				}
				case 2://Вправо
				{
					c[q].normals[h].set(sin(fi) / cos(fi), 1 / cos(fi), -sin(fi));
					qq2 = q2 == grainCount - 1 ? 0 : q2 + 1;
					y = 4;
					break;
				}
				case 3://На нас
				{
					c[q].normals[h].set(1 / cos(fi), -sin(fi), sin(fi) / cos(fi));
					qq1 = q1 == grainCount - 1 ? 0 : q1 + 1;
					y = 1;
					break;
				}
				case 4://Влево
				{
					c[q].normals[h].set(sin(fi), -1 / cos(fi), -sin(fi) / cos(fi));
					qq2 = q2 == 0 ? grainCount - 1 : q2 - 1;
					y = 2;
					break;
				}
				case 5://Вниз
				{
					c[q].normals[h].set(sin(fi) / cos(fi), sin(fi), -1 / cos(fi));
					qq3 = q3 == 0 ? grainCount - 1 : q3 - 1;
					y = 0;
					break;
				}
				//Далее идут уже необязательные соседи
				/**************           Рёбра куба          ***************************/
				case 6://Лево от нас
				{
					c[q].normals[h].set(-cos(fi + PI_2) * cos(fi)*cos(PI_2), -cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
					qq1 = q1 == 0 ? grainCount - 1 : q1 - 1;
					qq2 = q2 == 0 ? grainCount - 1 : q2 - 1;
					y = 9;
					break;
				}
				case 7://Лево на нас
				{
					c[q].normals[h].set(cos(fi + PI_2) * cos(fi)*cos(PI_2), -cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
					qq1 = q1 == grainCount - 1 ? 0 : q1 + 1;
					qq2 = q2 == 0 ? grainCount - 1 : q2 - 1;
					y = 8;
					break;
				}
				case 8://право от нас
				{
					c[q].normals[h].set(-cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
					qq1 = q1 == 0 ? grainCount - 1 : q1 - 1;
					qq2 = q2 == grainCount - 1 ? 0 : q2 + 1;
					y = 7;
					break;
				}
				case 9://Право на нас
				{
					c[q].normals[h].set(cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
					qq2 = q2 == grainCount - 1 ? 0 : q2 + 1;
					qq1 = q1 == grainCount - 1 ? 0 : q1 + 1;
					y = 6;
					break;
				}
				case 10://Верх лево
				{
					c[q].normals[h].set(cos(fi + PI_2) * sin(fi), -cos(fi + PI_2) * cos(fi), sin(fi + PI_2)*cos(fi));
					qq3 = q3 == grainCount - 1 ? 0 : q3 + 1;
					qq2 = q2 == 0 ? grainCount - 1 : q2 - 1;
					y = 15;
					break;
				}
				case 11://Верх право
				{
					c[q].normals[h].set(cos(fi + PI_2) * sin(fi), cos(fi + PI_2) * cos(fi), sin(fi + PI_2)*cos(fi));
					qq3 = q3 == grainCount - 1 ? 0 : q3 + 1;
					qq2 = q2 == grainCount - 1 ? 0 : q2 + 1;
					y = 14;
					break;
				}
				case 12://Верх на нас
				{
					c[q].normals[h].set(cos(fi + PI_2) * cos(fi), cos(fi + PI_2) * sin(fi), sin(fi + PI_2)*cos(fi));
					qq3 = q3 == grainCount - 1 ? 0 : q3 + 1;
					qq1 = q1 == grainCount - 1 ? 0 : q1 + 1;
					y = 17;
					break;
				}
				case 13://Верх от нас
				{
					c[q].normals[h].set(-cos(fi + PI_2) * cos(fi), cos(fi + PI_2) * sin(fi), sin(fi + PI_2)*cos(fi));
					qq3 = q3 == grainCount - 1 ? 0 : q3 + 1;
					qq1 = q1 == 0 ? grainCount - 1 : q1 - 1;
					y = 16;
					break;
				}
				case 14://Низ лево
				{
					c[q].normals[h].set(-cos(fi + PI_2) * sin(fi), -cos(fi + PI_2) * cos(fi), -sin(fi + PI_2)*cos(fi));
					qq2 = q2 == 0 ? grainCount - 1 : q2 - 1;
					qq3 = q3 == 0 ? grainCount - 1 : q3 - 1;
					y = 11;
					break;
				}
				case 15://Низ право
				{
					c[q].normals[h].set(-cos(fi + PI_2) * sin(fi), cos(fi + PI_2) * cos(fi), -sin(fi + PI_2)*cos(fi));
					qq2 = q2 == grainCount - 1 ? 0 : q2 + 1;
					qq3 = q3 == 0 ? grainCount - 1 : q3 - 1;
					y = 10;
					break;
				}
				case 16://Низ на нас
				{
					c[q].normals[h].set(cos(fi + PI_2) * cos(fi), -cos(fi + PI_2) * sin(fi), -sin(fi + PI_2)*cos(fi));
					qq1 = q1 == grainCount - 1 ? 0 : q1 + 1;
					qq3 = q3 == 0 ? grainCount - 1 : q3 - 1;
					y = 13;
					break;
				}
				case 17://Низ от нас
				{
					c[q].normals[h].set(-cos(fi + PI_2) * cos(fi), -cos(fi + PI_2) * sin(fi), -sin(fi + PI_2)*cos(fi));
					qq1 = q1 == 0 ? grainCount - 1 : q1 - 1;
					qq3 = q3 == 0 ? grainCount - 1 : q3 - 1;
					y = 12;
					break;
				}
				/**************      Вершины     *****************/
				case 18://верх лево от нас
				{
					c[q].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
					qq1 = q1 == 0 ? grainCount - 1 : q1 - 1;
					qq2 = q2 == 0 ? grainCount - 1 : q2 - 1;
					qq3 = q3 == grainCount - 1 ? 0 : q3 + 1;
					y = 25;
					break;
				}
				case 19://верх лево на нас
				{
					c[q].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
					qq1 = q1 == grainCount - 1 ? 0 : q1 + 1;
					qq2 = q2 == 0 ? grainCount - 1 : q2 - 1;
					qq3 = q3 == grainCount - 1 ? 0 : q3 + 1;
					y = 24;
					break;
				}
				case 20://верх право от нас
				{
					c[q].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
					qq1 = q1 == 0 ? grainCount - 1 : q1 - 1;
					qq2 = q2 == grainCount - 1 ? 0 : q2 + 1;
					qq3 = q3 == grainCount - 1 ? 0 : q3 + 1;
					y = 23;
					break;
				}
				case 21://верх право на нас
				{
					c[q].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
					qq1 = q1 == grainCount - 1 ? 0 : q1 + 1;
					qq2 = q2 == grainCount - 1 ? 0 : q2 + 1;
					qq3 = q3 == grainCount - 1 ? 0 : q3 + 1;
					y = 22;
					break;
				}
				case 22://низ лево от нас
				{
					c[q].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
					qq1 = q1 == 0 ? grainCount - 1 : q1 - 1;
					qq2 = q2 == 0 ? grainCount - 1 : q2 - 1;
					qq3 = q3 == 0 ? grainCount - 1 : q3 - 1;
					y = 21;
					break;
				}
				case 23://низ лево на нас
				{
					c[q].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
					qq1 = q1 == grainCount - 1 ? 0 : q1 + 1;
					qq2 = q2 == 0 ? grainCount - 1 : q2 - 1;
					qq3 = q3 == 0 ? grainCount - 1 : q3 - 1;
					y = 20;
					break;
				}
				case 24://низ право от нас
				{
					c[q].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));

					qq1 = q1 == 0 ? grainCount - 1 : q1 - 1;
					qq2 = q2 == grainCount - 1 ? 0 : q2 + 1;
					qq3 = q3 == 0 ? grainCount - 1 : q3 - 1;
					y = 19;
					break;
				}
				case 25://низ право на нас
				{
					c[q].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
					qq1 = q1 == grainCount - 1 ? 0 : q1 + 1;
					qq2 = q2 == grainCount - 1 ? 0 : q2 + 1;
					qq3 = q3 == 0 ? grainCount - 1 : q3 - 1;
					y = 18;
					break;
				}
				}
				int index = get1DPos(qq1, qq2, qq3);
				c[q].surrounds[h] = &c[index];					// Устанавливаем ссылки на соседние зерна
				c[index].surrounds[y] = &c[q];					// И в обратную сторону
				c[q].normals[h].normalize();

				for (int i = 0; i < DIM; i++)
				{
					c[index].normals[y].c[i] = -c[q].normals[h].c[i];// Сохранение нормалей граничащих фасеток
				}

				if (h < 6) c[q].contact[h] = 1;					//Контакт на грани октаэдра
				else if (h < 14) c[q].contact[h] = 3;			//Контакт на вершине октаэдра
				else c[q].contact[h] = 2;						//Контакт на ребре
			}
			if (prms::grainSurroundCount > 6)					//Уменьшение объёма зерна из-за отсечений
			{
				// Аппроксимация объема, уменьшающегося из-за отсечения
				// вершин и ребер. Не претендует на точность

				double a = c[q].size * 0.1;						//Длина срезанной части вдоль ребра
				double volEdge = a * a * c[q].size / 2.0;		//Объём, срезанный рёбрами
				double volVertex = a * a * a / SQRT3;			//Объём, срезанный вершинами
				int cut_edge = 0;								//Кол-во срезанных рёбер
				int cut_vertex = 0;								//Кол-во срезанных вершин
				for (int h = 6; h < prms::grainSurroundCount; h++)
				{
					if (c[q].contact[h] != 0)
					{
						if (h < 14) cut_vertex++;
						else cut_edge++;
					}
				}
				c[q].volume -= (cut_edge*volEdge + cut_vertex * volVertex);//Вычитание срезанного объема
			}


		}
	}
}