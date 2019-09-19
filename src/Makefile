A_TARGET=snap

# -xAVX -pg g 
# C compiler command
CC=mpiicc
CFLAGS=-mkl -O3 -std=c11 -xAVX -fopenmp -parallel -DUSEMKL -DUSEVML -pg -g # C compiler flags
LDFLAGS=-L/usr/local/cuda/lib64 # Linker flags, e.g. -L<lib dir> if you have libraries in a nonstandard directory <lib dir>
LIBS=-lcudart # Libraries to pass to the linker, e.g. -l<library>

# C++
CXX=icpc # C++ compiler command
CXXFLAGS=-O3 -std=c++11 -qopenmp -pg -g	# C++ compiler flags
CXXCPP= # C++ preprocessor flags

# GPU
GPUCC=nvcc
GPUFLAGS=-dc -arch=sm_60 -dlink -m64 -rdc=true


A_FILES := $(wildcard *.cpp) $(wildcard *.c) $(wildcard *.cu)
A_OBJS := $(sort $(patsubst %.c, %.o, $(patsubst %.cu, %.o, $(patsubst %.cpp, %.o, $(A_FILES)))))


SRCS = snap.h snap_main.c utils.c version.c plib.c geom.c \
       sn.c data.c control.c input.c setup.c dealloc.c \
       translv.c solvar.c outer.c expxs.c inner.c \
       octsweep.c dim1_sweep.c dim3_sweep.c output.c time.c mms.c


all: $(A_TARGET) 			# Build all targets

vars:					# Show Makefile variable values
	@echo "       CC "$(CC)
	@echo " A_TARGET "$(A_TARGET)
	@echo "   CFLAGS "$(CFLAGS)
	@echo "  LDFLAGS "$(LDFLAGS)
	@echo "     LIBS "$(LIBS)@egrep "[a-z]\w+:" Makefile | sort
	@echo "   A_OBJS "$(A_OBJS)

clean:					# Remove make build object and target files
	-rm -f *.o $(A_TARGET)

cleanout: clean
	-rm -f $(OUTPUT_FILE)

.PHONY:	clean vars help rebuild

# Clean and build target
rebuild: clean $(A_TARGET)

# Build target and its required object files
$(A_TARGET): $(A_OBJS)
	@echo "Compiling and linking main program"
	$(GPUCC) -arch=sm_60 -m64 -std=c++11 -ccbin=$(CC) -Xcompiler "-mkl" $(LDFLAGS) $(LIBS) $^ -o $@

%.o: %.cu
	$(GPUCC) $(GPUFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CXXCPP) -c $< -o $@

%: %.cpp
	$(CXX) $(CXXFLAGS) $(CXXCPP) -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

%: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $<

# Count lines of code
count:
	rm -f Lines
	for file in $(SRCS); do ./LineCount $$file Lines; done
	gawk -f ./LineReport < Lines >> Lines
	cat Lines
