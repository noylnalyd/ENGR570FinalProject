temp:	Simulator.hpp Simulator.cpp washer.hpp BodyModel.hpp BodyModel.cpp element.cpp element.hpp PseudoBlockMatrix.hpp sector.hpp SimModel.hpp
	gcc -o temp.exe Simulator.cpp -I. -g
runtemp:
	./temp.exe