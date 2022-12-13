temp:	washer.hpp sector.hpp element.hpp BodyModel.hpp PseudoBlockMatrix.hpp SimModel.hpp Simulator.hpp element.cpp BodyModel.cpp Simulator.cpp
	g++ -o temp.exe Simulator.cpp -g -std=c++11

test:	washer.hpp sector.hpp element.hpp BodyModel.hpp PseudoBlockMatrix.hpp SimModel.hpp Simulator.hpp element.cpp BodyModel.cpp Simulator.cpp
	g++ test.cpp -o test.exe -g -std=c++11
runtemp:
	./temp.exe
runtest:
	./test.exe