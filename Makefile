FC = gfortran
FCFLAG=

# for gfortran
ifeq (${FC},gfortran)
MODFLAG= -J
FCFLAG+= -fopenmp -ffree-line-length-512 -ffpe-trap=invalid,zero,overflow -fbacktrace #-fmax-stack-var-size=32768 
endif

# for ifort
ifeq (${FC},ifort)
MODFLAG= -module 
FCFLAG+= -qopenmp #-heap-arrays -traceback -fpe0 -CB #-pg -check all 
endif

OBJ_DIR = $(HORMONE_DIR)/obj
SRC_DIR = $(HORMONE_DIR)/src
BIN_DIR = $(HORMONE_DIR)/bin

TARGET = hormone

# compile modules that are being depended on first
MOD = modules.f90 recombination.f90 eos.f90 fluxlimiter.f90 hlldflux.f90 dirichlet.f90 particles.f90 cooling.f90
MODF= $(patsubst %,$(SRC_DIR)/%,$(MOD))
TARGETF= $(patsubst %,$(SRC_DIR)/%.f90,$(TARGET))
SRC = $(MODF) \
      $(filter-out $(MODF) $(TARGETF), $(sort $(wildcard $(SRC_DIR)/*.f90))) \
      $(TARGETF)
OBJ = $(patsubst $(SRC_DIR)/%.f90,%.o,$(SRC))
OBJ_F = $(patsubst %.o,$(OBJ_DIR)/%.o,$(OBJ))

all: $(TARGET)

$(TARGET): $(OBJ_F)
	$(FC) $(FCFLAG) -o $(BIN_DIR)/$@ $(OBJ_F)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FCFLAG) -c $< -o $@ $(MODFLAG)$(OBJ_DIR)

install: $(TARGET)
	install -s $(TARGET) $(BIN_DIR)

.PHONY: clean
clean:
	rm -f $(OBJ_DIR)/*.o *~ $(SRC_DIR)/*~ $(OBJ_DIR)/*.mod $(BIN_DIR)/*