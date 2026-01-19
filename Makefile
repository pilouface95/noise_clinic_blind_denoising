# Makefile created with Create Makefile.sh 18/10/2013
# (c) Miguel Colom
# (c) Marc Lebrun

COPT      = -O99
CFLAGS    = $(COPT) -fopenmp -ggdb
CSTRICT   = -Wall -Wextra -ansi
CXXOPT    = -O99
CXXFLAGS  = $(CXXOPT) -fopenmp -ggdb
CXXSTRICT = -Wall -Wextra -ansi
LDFLAGS   = -lpng -lgomp -lfftw3f -lfftw3f_threads
EXEC      = MSD

default: $(EXEC)
all: $(EXEC)

# ------- C++ files -------
./NoiseEstimate/RunNoiseEstimate.o: ./NoiseEstimate/RunNoiseEstimate.cpp ./NoiseEstimate/RunNoiseEstimate.h
	$(CXX) $(CXXFLAGS) -c ./NoiseEstimate/RunNoiseEstimate.cpp -o ./NoiseEstimate/RunNoiseEstimate.o
	
./Utilities/CImage.o: ./Utilities/CImage.cpp ./Utilities/CImage.h
	$(CC) $(CFLAGS) -c ./Utilities/CImage.cpp -o ./Utilities/CImage.o

./NoiseEstimate/CurveFiltering.o: ./NoiseEstimate/CurveFiltering.cpp ./NoiseEstimate/CurveFiltering.h
	$(CXX) $(CXXFLAGS) -c ./NoiseEstimate/CurveFiltering.cpp -o ./NoiseEstimate/CurveFiltering.o

./NoiseEstimate/CHistogram.o: ./NoiseEstimate/CHistogram.cpp
	$(CXX) $(CXXFLAGS) -c ./NoiseEstimate/CHistogram.cpp -o ./NoiseEstimate/CHistogram.o


./NoiseEstimate/	.o: ./NoiseEstimate/operations.cpp
	$(CXX) $(CXXFLAGS) -c ./NoiseEstimate/operations.cpp -o ./NoiseEstimate/operations.o

./main.o: ./main.cpp
	$(CXX) $(CXXFLAGS) -c ./main.cpp -o ./main.o

./LibNlBayes/NlBayes.o: ./LibNlBayes/NlBayes.cpp ./LibNlBayes/NlBayes.h
	$(CXX) $(CXXFLAGS) -c ./LibNlBayes/NlBayes.cpp -o ./LibNlBayes/NlBayes.o

./LibNlBayes/LibMatrix.o: ./LibNlBayes/LibMatrix.cpp ./LibNlBayes/LibMatrix.h
	$(CXX) $(CXXFLAGS) -c ./LibNlBayes/LibMatrix.cpp -o ./LibNlBayes/LibMatrix.o

./Utilities/Utilities.o: ./Utilities/Utilities.cpp ./Utilities/Utilities.h
	$(CXX) $(CXXFLAGS) -c ./Utilities/Utilities.cpp -o ./Utilities/Utilities.o

./Utilities/LibImages.o: ./Utilities/LibImages.cpp ./Utilities/LibImages.h
	$(CXX) $(CXXFLAGS) -c ./Utilities/LibImages.cpp -o ./Utilities/LibImages.o

./MultiScalesDenoising/MultiScalesDenoising.o: ./MultiScalesDenoising/MultiScalesDenoising.cpp ./MultiScalesDenoising/MultiScalesDenoising.h
	$(CXX) $(CXXFLAGS) -c ./MultiScalesDenoising/MultiScalesDenoising.cpp -o ./MultiScalesDenoising/MultiScalesDenoising.o

# ------- Main -------
$(EXEC):  ./Utilities/CImage.o ./NoiseEstimate/RunNoiseEstimate.o ./NoiseEstimate/CurveFiltering.o ./main.o ./LibNlBayes/NlBayes.o ./LibNlBayes/LibMatrix.o ./Utilities/Utilities.o ./Utilities/LibImages.o ./MultiScalesDenoising/MultiScalesDenoising.o
	$(CXX)  ./Utilities/CImage.o ./NoiseEstimate/RunNoiseEstimate.o ./NoiseEstimate/CurveFiltering.o ./main.o ./LibNlBayes/NlBayes.o ./LibNlBayes/LibMatrix.o ./Utilities/Utilities.o ./Utilities/LibImages.o ./MultiScalesDenoising/MultiScalesDenoising.o $(LDFLAGS) -o $(EXEC)

lint: 
	$(MAKE) CFLAGS="$(CFLAGS) $(CSTRICT)" CXXFLAGS="$(CXXFLAGS) $(CXXSTRICT)"

clean: 
	rm -f *.o

distclean: clean
	rm -f $(EXEC)

