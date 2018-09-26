// This is an open source non-commercial project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include "stdafx.h"

#include "Params.h"
#include "GrainStructure.h"

namespace model
{

	double rnd() { return double(rand()) / RAND_MAX; }

	GrainStructure::GrainStructure(Polycrystall* poly, bool periodic)
	{
		this->polycrystall = poly;
		this->periodic = periodic;
	}


	voro::container* GrainStructure::makeContainer()
	{
		// Количество вычислительных блоков для разбиения контейнера
		const int nx = 6, ny = 6, nz = 6;
		// Задание контейнера с описанной выше геометрией
		// Периодическими ГУ по всем направлениям
		// И выделяющий память под 8 частиц в ккаждом блоке
		return new voro::container(xMin, xMax, yMin, yMax, zMin, zMax, nx, ny, nz,
			periodic, periodic, periodic, 8);
	}

	// Укладка зерен в виде многогранников Вороного. Разбиение моделируемого объема среды
	// на ячейки без пустот и наложений. Используется пакет voro++
	void GrainStructure::makeVoronoiStructure()
	{
		voro::container *con = makeContainer();
		// Количество частиц - центров кристаллизации
		const int particles = polycrystall->totalGrainCount;
		
		// Добавление случайно сгенерированных частиц в контейнер
		double x, y, z;
		for (int i = 0; i < particles; i++) {
			x = xMin + rnd()*(xMax - xMin);
			y = yMin + rnd()*(yMax - yMin);
			z = zMin + rnd()*(zMax - zMin);
			con->put(i, x, y, z);
			posMap.insert(std::pair<int, Vector>(i, Vector(x, y, z)));
			//centerPos.push_back(Vector(x, y, z));
			//ids.push_back(i);
		}
		updateStructure(con);
	}

	void GrainStructure::updateContainer()
	{
		voro::container *con = makeContainer();

		double x, y, z;

		for (std::pair<int, Vector> pair : posMap)
		{
			Vector v = pair.second;
			con->put(pair.first, v.c[0], v.c[1], v.c[2]);
		}
	/*	for (int i = 0; i < posMap.size(); i++)
		{
			x = posMap[i].c[0];
			y = posMap[i].c[1];
			z = posMap[i].c[2];
			int index = ids[i];
			con->put(index, x, y, z);
		}*/

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
		std::map<int, Vector>::iterator it;
		it = posMap.begin();
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

		//	Grain* grain = polycrystall->findGrainById(ids[q]);
			int pos = polycrystall->findGrainPosById(it->first);
			// Выделение памяти для всех топологических параметров структуры
		//	grain->clearTopology();
			polycrystall->c[pos].neighbors = std::vector<Grain*>(neighborCount);
			polycrystall->c[pos].normals = std::vector<Vector>(neighborCount);
			polycrystall->c[pos].areas = std::vector<double>(neighborCount);

			for (int i = 0; i < neighborCount; i++)
			{
				int curr = i * 3;
				Vector normal(normals[curr], normals[curr + 1], normals[curr + 2]);
				// Нормали уже отнормированы
				polycrystall->c[pos].normals[i] = normal;
				polycrystall->c[pos].neighbors[i] = &polycrystall->c[neighbors[i]];
				polycrystall->c[pos].areas[i] = areas[i];
			}
			// Начальный объем каждого зерна
			polycrystall->c[pos].volume = pow(polycrystall->c[q].size, 3) * cell.volume();
			polycrystall->c[pos].position = q;
			it++;
			q++;
		} while (loop.inc());

		lastId = q;
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