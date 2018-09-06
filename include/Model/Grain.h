#pragma once
#ifndef __FRAGMENT_H 
#define __FRAGMENT_H

#include "vector"

#include "Tensor.h"
#include "Tensor4.h"
#include "SlipSystem.h"

namespace model
{

	class Grain
	{
	public:
		Tensor d;							// Тензор деформации скорости
		Tensor w;							// Тензор вихря
		Tensor d_in;						// Тензор неупругой составляющей деформации
		Tensor o;							// Ориентационный тензор
		Tensor om;							// Тензор спина решётки
		Tensor sgm;							// Тензор напряжений
		Tensor dsgm;						// Тензор скоростей напряжений
		Tensor e;							// Тензор деформаций
		Tensor4 p;							// Тензор упругих свойств
		
		int ssCount;						// Кол-во систем скольжения
		SlipSystem *ss;						// Системы скольжения

		int materialType;					// Тип материала
		double size;						// Линейный размер фрагмента
		double volume;						// Объём фрагмента с учётом срезов
		double stress;						// Интенсивность напряжений
		double strain;						// Интенсивность деформаций
		
		std::vector<Grain*> neighbors;		// Ссылки на граничащие фрагменты
		std::vector<Vector> normals;		// Вектора нормали к граничащим фрагментам
		std::vector<double> areas;			// Площади фасеток
		Vector moment;								

		int position;						// Порядковый номер данного элемента в поликристалле

		bool isRotate;						// Вращается ли решётка на текущем шаге
		double rotationSpeed;				// Скорость вращения решётки
		double rotationTotalAngle;			// Накопленный угол поворота решётки
		double rotationEnergy;				// Энергия ротаций фрагмента
		
		double rotationParamMc;				// Начальный критический момент
		double rotationParamH;				// Коэффициент при необратимой составляющей поворота
		double rotationParamA;				// Коэффициент при обратимой составляющей поворота
		double rotationParamL;				// Коэффициент при скорости поверхностных моментов

		// Задание материальных параметров
		void setMaterialParams(int);
		// Задание начальной ориентации по углам Эйлера
		void orientate(double, double, double);
		// Задание начальной ориентации по оси и углу
		void orientateAxis(double, Vector);
		// Задание начальной ориентации кватернионом
		void orientateQuater(double, double, double, double);
		// Вычисление НДС фрагмента
		void stressStrainCalc();
		// Возвращает меру разориентации с выбранным соседом
		double disorientMeasure(int);

		Grain();
		~Grain();

	private:

	};

}

#endif