// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"

#include "Params.h"
#include "GrainStructure.h"

namespace model
{
		
	double rnd() { return double(rand()) / RAND_MAX; }

	GrainStructure::GrainStructure(Polycrystall* poly)
	{
		this->polycrystall = poly;
	}

	// Укладка зерен в виде многогранников Вороного. Разбиение моделируемого объема среды
	// на ячейки без пустот и наложений. Используется пакет voro++
	void GrainStructure::makeVoronoiStructure()
	{
		// Геометрия контейнера для расчетной области
		const double x_min = 0, x_max = 1;
		const double y_min = 0, y_max = 1;
		const double z_min = 0, z_max = 1;
		// Количество вычислительных блоков для разбиения контейнера
		const int n_x = 6, n_y = 6, n_z = 6;
		// Количество частиц - центров кристаллизации
		const int particles = polycrystall->totalGrainCount;
		// Задание контейнера с описанной выше геометрией
		// Периодическими ГУ по всем направлениям
		// И выделяющий память под 8 частиц в ккаждом блоке
		voro::container *con = new voro::container(x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z,
			true, true, true, 8);
		centerPos.clear();
		// Добавление случайно сгенерированных частиц в контейнер
		double x, y, z;
		for (int i = 0; i < particles; i++) {
			x = x_min + rnd()*(x_max - x_min);
			y = y_min + rnd()*(y_max - y_min);
			z = z_min + rnd()*(z_max - z_min);
			con->put(i, x, y, z);
			centerPos.push_back(Vector(x, y, z));
		}
		updateStructure(con);
	}

	void GrainStructure::updateStructure(voro::container* con)
	{
		// Ячейка, к которой будем обращаться в цикле
		voro::voronoicell_neighbor cell;
		// Вспомогательный объект, обходящий все элементы контейнера
		voro::c_loop_all loop(*con);
		// Основной цикл обхода контейнера
		int q = 0;
		if (loop.start()) do if (con->compute_cell(cell, loop))
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
			polycrystall->c[q].neighbors = std::vector<Grain*>(neighborCount);
			polycrystall->c[q].normals = std::vector<Vector>(neighborCount);
			polycrystall->c[q].areas = std::vector<double>(neighborCount);

			for (int i = 0; i < neighborCount; i++)
			{
				int curr = i * 3;
				Vector normal(normals[curr], normals[curr + 1], normals[curr + 2]);
				// Нормали уже отнормированы
				polycrystall->c[q].normals[i] = normal;
				polycrystall->c[q].neighbors[i] = &polycrystall->c[neighbors[i]];
				polycrystall->c[q].areas[i] = areas[i];
			}
			// Начальный объем каждого зерна
			polycrystall->c[q].volume = pow(polycrystall->c[q].size, 3) * cell.volume();
			polycrystall->c[q].position = q;
			q++;
		} while (loop.inc());
		int a = 9;
	}

	void printStructureInfo(voro::container *con)
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
		(*con).print_custom("%i %s %t %l %n %v", "structure.txt");
	}
}