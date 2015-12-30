NAME = $(shell uname -n)
SHELL_V = $(shell uname -s)_Kernel_$(shell uname -v | awk '{print $$4}' | sed s/":"//)
ARCH = $(shell uname -p)
BITS = $(shell uname -m)
OS = $(shell uname -s)

USERMACROS =
CXXFLAGS =
DFLAGS =
LFLAGS =
INSTALL =
EXE_BITS =
EXE_OS =
#HAS_FFT_LIB =$(strip $(wildcard ../lib/libfftw3.a))
#HAS_MPI =$(strip $(wildcard /usr/local/lib/libmpi.so))
#HAS_BLAS_LIB =$(strip $(wildcard ../lib/libblas.so.3))
#HAS_CLAPACK_LIB =$(strip $(wildcard ../lib/liblapack.a))
#HAS_MINPACK_LIB =$(strip $(wildcard ../lib/libcminpack.a))
#HAS_CUBATURE_LIB =$(strip $(wildcard ../lib/libcubature.a))

SRC := union.cpp stopwatch.cpp
OBJS := union.o stopwatch.o

OUTLIB := libunionfind.a

#====================================
#  makefile options
#  
# compiler
ifeq ($(COMP),)
   CXX = g++
else
   CXX = $(COMP)
endif

ifeq ($(OS),Linux)
   EXE_OS = nix
   USERMACROS += -DOS=\"unix\"
endif

# make shared lib
ifeq ($(MAKE_LIB),1)
   USERMACROS += -fPIC
endif

#mode
ifeq ($(MODE),)
   MODE = debug
   USERMACROS += -DDEBUG
   CXXFLAGS += -Wall
   CXXFLAGS += -g
   CXXFLAGS += -ggdb
else
   MODE = release
   CXXFLAGS += -g
   CXXFLAGS += -O3
   USERMACROS += -DNDEBUG	
endif

#bits
ifeq ($(BITS),x86_64)
   EXE_BITS = 64
   CXXFLAGS += -m64
   USERMACROS += -DBIT_64
else
   EXE_BITS = 32
   CXXFLAGS += -m32
   USERMACROS += -DBIT_32
endif

# threads
ifeq ($(THREADED),)
   USERMACROS += -DTHREADED
   CXXFLAGS += -pthread
   LFLAGS += -lpthread 
endif

# cubature integration
ifeq ($(USE_HCUBATURE),true)
  USERMACROS += -DHCUBA
  LFLAGS += -lm
  INCS +=-I../external/cubature
endif

ifeq ($(HAS_MINPACK_LIB),)
else
   USERMACROS += -DHAS_MINPACK
   LFLAGS +=-L../lib -lcminpack
   INCS +=-I../external/cminpack
endif

# check for fftw
ifeq ($(HAS_FFT_LIB),)
   BUILD_FFT = fftw
else
   USERMACROS += -DHAS_FFTW
   LFLAGS +=-L../lib -lfftw3
endif

# check for BLAS
ifeq ($(HAS_BLAS_LIB),)
   BUILD_BLAS = blas
else
   USERMACROS += -DHAS_BLAS
   LFLAGS +=-L/usr/lib -lblas
   INCS +=-I../external/cblas/include
endif

# check for lapack
ifeq ($(HAS_CLAPACK_LIB),)
   BUILD_CLAPACK = lapack
else
   USERMACROS += -DHAS_LAPACK
   LFLAGS +=-L../lib -llapack
   INCS +=-I../external/lapack/INCLUDE
endif


# executable 
EXE = unionfind-$(EXE_OS)-$(EXE_BITS)

# date..
USERMACROS += -DBUILD_DATE="\"$$(date)\""
USERMACROS += -DVERSION="\"ver-1.0-$$(date)\""
# check for mpi
ifeq ($(HAS_MPI),)
   NO_MPI = true
else
   USERMACROS += -DHAS_MPI
#   LFLAGS += -lopenmpi
   CXX = mpicxx 
endif

# includes
CXXFLAGS +=$(INCS)

# collect object files here
OBJ := $(patsubst %.cpp, %.o, $(filter %.cpp,$(SRC)))

.PHONY:all

all: clean information $(OUTLIB) copy $(BUILD_CONTROL) clean


information:
	@echo ""
	@echo "+------------------------------------+"
	@echo "|       -- Build info -- "
	@echo "|..NAME	     = "$(NAME)
	@echo "|..SHELL     = "$(SHELL_V)
	@echo "|..MODE      = "$(MODE)
	@echo "|..ARCH      = "$(ARCH) $(EXE_BITS)"-bit"
	@echo "|..OS        = "$(EXE_OS)
	@echo "|..MPI       = "$(HAS_MPI)
	@echo "|..EXE       = "$(EXE)
	@echo "+------------------------------------+"
	@echo ""

#linking the program
link: $(OBJ)
	$(CXX) -o $(EXE) $(OBJS) $(LFLAGS)	

#$(CXX) $(USERMACROS) $(CXXFLAGS) -c $< -o $@

%.o:%.cpp
	$(CXX) $(USERMACROS) $(CXXFLAGS) -c $< $(LFLAGS)

$(OUTLIB): $(OBJS)
	ar ru $@ $^
	ranlib $@

copy:
	mv $(OUTLIB) ../lib

TAGS:   $(SRC)
	etags $(SRC)

.PHONY:clean
clean:
	@echo ""
	@echo "..Clean source directory.."
	@echo ""
	find . -name "*.o" | xargs rm -vf
	rm -vf *~ *#
