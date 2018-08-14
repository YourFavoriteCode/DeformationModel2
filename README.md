# DeformationModel
Статистическая двухуровневая математическая модель неупругого деформирования представительного объема поликристалла

Модель позволяет определять отклик материала при различных видах нагружения, строить кривые напряженно-деформированного состояния (НДС), строить прямые полюсные фигуры (ППФ), проводить анализ любых внутренних переменных.
В основе модели лежит гипотеза Фойгта. Элементом верхнего масштабного уровня (макро-) является представительный объем поликристалла. На нижнем (мезо-) уровне рассматриваются зерна либо их фрагменты.
Определяющим соотношением является закон Гука в скоростной форме. Для определения скоростей сдвигов дислокаций используется закон Хатчинсона. 

# Особенности
Реализованы механизмы деформационного упрочнения (базовое изотропное и ориентированное зернограничное слагаемые).
Реализованы модели ротаций кристаллических решеток по модели Тейлора, а также по модели, учитывающей несовместность пластических деформаций.
Возможна работа с ГЦК и ОЦК материалами, а также со смешанными материалами.

# Сборка
Проект полностью написан на  C++.
Требуется Visual Studio (2012 и выше). Допускается распараллеливание вычислений с помощью технологии OpenMP (нужен компилятор с ее поддержкой, например Intel C++ compiler).
Зависимости: Eigen (https://github.com/eigenteam/eigen-git-mirror) для решения СЛАУ и TinyXML (https://github.com/leethomason/tinyxml2) для удобного конфигурирования модели.

# Конфигурирование
Настройка всех параметров модели осуществляется путем их задания в XML-файле, лежащем в папке с исполняемым файлом модели, и вызова исполняемого файла с единственным параметром командной строки - названием XML-файла.
