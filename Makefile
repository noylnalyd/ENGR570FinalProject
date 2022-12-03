temp:	washer.hpp sector.hpp element.hpp BodyModel.hpp PseudoBlockMatrix.hpp SimModel.hpp Simulator.hpp element.cpp BodyModel.cpp Simulator.cpp
	g++ -o temp.exe Simulator.cpp -g

test:	washer.hpp sector.hpp element.hpp BodyModel.hpp PseudoBlockMatrix.hpp SimModel.hpp Simulator.hpp element.cpp BodyModel.cpp Simulator.cpp
	g++ test.cpp -o test.exe -g
runtemp:
	./temp.exe