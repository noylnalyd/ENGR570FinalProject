test:	washer.hpp sector.hpp element.hpp BodyModel.hpp PseudoBlockMatrix.hpp SimModel.hpp Simulator.hpp element.cpp BodyModel.cpp Simulator.cpp
	g++ test.cpp -o test.exe -g -std=c++11
runtest:
	./test.exe