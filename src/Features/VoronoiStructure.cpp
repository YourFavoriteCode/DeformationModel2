// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"

#include "Params.h"
#include "GrainStructure.h"

#include "voro++.hh"

using namespace voro;

namespace model
{
		
	double rnd() { return double(rand()) / RAND_MAX; }

	// Укладка зерен в виде многогранников Вороного. Разбиение моделируемого объема среды
	// на ячейки без пустот и наложений. Используется пакет voro++
	void GrainStructure::makeVoronoiStructure(Polycrystall* poly)
	{
		// Геометрия контейнера для расчетной области
		const double x_min = 0, x_max = 1;
		const double y_min = 0, y_max = 1;
		const double z_min = 0, z_max = 1;
		// Количество вычислительных блоков для разбиения контейнера
		const int n_x = 6, n_y = 6, n_z = 6;
		// Количество частиц - центров кристаллизации
		const int particles = poly->totalGrainCount;
		// Задание контейнера с описанной выше геометрией
		// Периодическими ГУ по всем направлениям
		// И выделяющий память под 8 частиц в ккаждом блоке
		container con(x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z,
			true, true, true, 8);
		centerPos.clear();
		// Добавление случайно сгенерированных частиц в контейнер
		double x, y, z;
		for (int i = 0; i < particles; i++) {
			x = x_min + rnd()*(x_max - x_min);
			y = y_min + rnd()*(y_max - y_min);
			z = z_min + rnd()*(z_max - z_min);
			con.put(i, x, y, z);
			centerPos.push_back(Vector(x, y, z));
		}
		// Ячейка, к которой будем обращаться в цикле
		voronoicell_neighbor cell;
		// Вспомогательный объект, обходящий все элементы контейнера
		c_loop_all loop(con);
		// Основной цикл обхода контейнера
		int q = 0;
		if (loop.start()) do if (con.compute_cell(cell, loop))
		{
			std::vector<double> normals;	// Нормали к каждой фасетке
			std::vector<double> areas;		// Плаощади фасеток
			std::vector<int> neighbors;		// ID соседнего к каждой фасетке зерна
			std::vector<int> vertices;		// Порядок вершин многогранника
			cell.neighbors(neighbors);
			cell.normals(normals);
			cell.face_areas(areas);
			cell.vertex_orders(vertices);

			int neighborCount = neighbors.size();
			int verticesCount = vertices.size();

			// Выделение памяти для всех топологических параметров структуры
			poly->c[q].neighbors = std::vector<Grain*>(neighborCount);
			poly->c[q].normals = std::vector<Vector>(neighborCount);
			poly->c[q].areas = std::vector<double>(neighborCount);

			for (int i = 0; i < neighborCount; i++)
			{
				int curr = i * 3;
				Vector normal(normals[curr], normals[curr + 1], normals[curr + 2]);
				// Нормали уже отнормированы
				poly->c[q].normals[i] = normal;
				poly->c[q].neighbors[i] = &poly->c[neighbors[i]];
				poly->c[q].areas[i] = areas[i];
			}
			// Начальный объем каждого зерна
			poly->c[q].volume = pow(poly->c[q].size, 3) * cell.volume();
			poly->c[q].position = q;
			q++;
		} while (loop.inc());
	}

	void printStructureInfo(container con)
	{
		
		// Легенда:
		// %i - ID ячейки
		// %q - координаты в формате %x %y %z
		// %w - кол-во вершин многогранника
		// %p - список вершин в формате %x %y %z относительно центра многогранника
		// %P - список вершин в формате %x %y %z в глобальной системе координат
		// %o - упорядоченный список вершин
		// %g - кол-во ребер ячейки
		// %E - сумма длин ребер
		// %e - список периметров каждой фасетки
		// %s - кол-во фасеток ячейки
		// %F - полная площадь поверхности ячейки
		// 
		//
		// 
		// %t - список вершин в скобках, образующие каждую грань в глобальной нумерации
		// %l - список нормалей ко всем фасеткам
		// %n - список соседних фасеток
		// %v - объем многогранника
		con.print_custom("%i %s %t %l %n %v", "structure.txt");
	}
}