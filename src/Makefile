A_TARGET=snap

# C compiler command
CC=mpiicc
# C compiler flags
CFLAGS=-O3 -std=c11 -xAVX -qopenmp -parallel -pg -g
# Linker flags, e.g. -L<lib dir> if you have libraries in a nonstandard directory <lib dir>
LDFLAGS=-L/usr/local/cuda/lib64
# Libraries to pass to the linker, e.g. -l<library>
LIBS=-lcudart -lcuda

# C++
# C++ compiler command
CXX=icpc
# C++ compiler flags
CXXFLAGS=-O3 -std=c++11 -qopenmp -pg -g
# C++ preprocessor flags
CXXCPP=

# GPU
GPUCC=nvcc
GPUFLAGS=-arch=compute_61 -code=sm_61 -m 64 -rdc=true -dc -dlink


A_FILES := $(wildcard *.cpp) $(wildcard *.c) $(wildcard *.cu)
A_OBJS := $(sort $(patsubst %.c, %.o, $(patsubst %.cu, %.o, $(patsubst %.cpp, %.o, $(A_FILES)))))


SRCS = snap.h snap_main.c utils.c version.c plib.c geom.c \
       sn.c data.c control.c input.c setup.c dealloc.c \
       translv.c solvar.c outer.c expxs.c inner.c \
       octsweep.c dim1_sweep.c output.c time.c mms.c

C_FILES := $(filter-out dim3_sweep.c, $(wildcard *.c))

# Build all targets
all: $(A_TARGET)

# Show Makefile variable values
vars:
	@echo "       CC "$(CC)
	@echo " A_TARGET "$(A_TARGET)
	@echo "   CFLAGS "$(CFLAGS)
	@echo "  LDFLAGS "$(LDFLAGS)
	@echo "     LIBS "$(LIBS)@egrep "[a-z]\w+:" Makefile | sort
	@echo "   A_OBJS "$(A_OBJS)
	@echo " C_FILES"$(C_FILES)

# Remove make build object and target files
clean:
	-rm -f *.o $(A_TARGET)

cleanout: clean
	-rm -f $(OUTPUT_FILE)

.PHONY:	clean vars help rebuild

# Clean and build target
rebuild: clean $(A_TARGET)

# Build target and its required object files
$(A_TARGET): $(A_OBJS)
	$(GPUCC) -arch=compute_61 -code=sm_61 -m 64 -std=c++11 -ccbin=$(CC) -Xcompiler "-liomp5" $(LDFLAGS) $(LIBS) $^ -o $@

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
