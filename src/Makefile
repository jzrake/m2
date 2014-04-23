# if there is no Makefile.in then use the template
# --------------------------------------------------
ifneq ($(strip $(MAKEFILE_IN)),)
# use value of MAKEFILE_IN if provided on the command line
else ifeq ($(shell test -e Makefile.in && echo 1), 1)
MAKEFILE_IN = $(shell pwd)/Makefile.in
else
MAKEFILE_IN = $(shell pwd)/Makefile.in.template
endif
include $(MAKEFILE_IN)


# git versioning
# --------------------------------------------------
GIT_SHA := $(shell git rev-parse HEAD | cut -c 1-10)
TMP_H := $(shell mktemp -u make.XXXXXX)


# m2 executables
# --------------------------------------------------
EXE = m2


# object code required for executables
# --------------------------------------------------
SRC1 = $(wildcard *.c)
SRC2 = $(wildcard *.cpp)
OBJ = $(SRC1:.c=.o) $(SRC2:.cpp=.o)
DEP = $(SRC1:.c=.dep) $(SRC2:.cpp=.dep)


LUA_L  = $(LUA_LIB)
ZULU_L = zulu/libzulu.a
GLUI_L = glui/libglui.a
INC = -Izulu -Ilua

LIB += $(LUA_L)
LIB += $(ZULU_L)

LIB_STATIC += $(LUA_L)
LIB_STATIC += $(ZULU_L)

ifeq ($(HAVE_GLUT), 1)
LIB += $(GLUT_L) $(GLUI_L)
INC += -Iglui
LIB_STATIC += $(GLUI_L)
endif

ifeq ($(HAVE_HDF5), 1)
LIB += -L$(HDF5_HOME)/lib -lhdf5
endif



# build rules
# --------------------------------------------------
default : $(EXE)

vis.o : vis.cpp m2-cfg.h $(MAKEFILE_IN)
	@$(CC) -MM $< $(INC) $(GLUT_I) > $(<:.cpp=.dep)
	$(CC) $(CFLAGS) $(INC) $< -c $(GLUT_I)

%.o : %.c m2-cfg.h $(MAKEFILE_IN)
	@$(CC) -MM $< $(INC) > $(<:.c=.dep)
	$(CC) $(CFLAGS) $(INC) $< -c

$(EXE) : $(OBJ) $(LIB_STATIC)
	$(CXX) $(CFLAGS) $(OBJ) $(LIB) -o $@

m2-cfg.h : .FORCE
	@$(RM) $(TMP_H)
	@echo "/* m2 config header file */" >> $(TMP_H)
	@echo "#define M2_HAVE_MPI $(HAVE_MPI)" >> $(TMP_H)
	@echo "#define M2_HAVE_GLUT $(HAVE_GLUT)" >> $(TMP_H)
	@echo "#define M2_GIT_SHA \"$(GIT_SHA)\"" >> $(TMP_H)
	@echo "#define M2_DEBUG_FILE $(DEBUG_FILE)" >> $(TMP_H)
	@echo "#define M2_INSTALL_PATH \"$(PWD)\"" >> $(TMP_H)
	@cmp -s $(TMP_H) $@ || (echo "[m2-cfg.h updated]"; cat $(TMP_H)>$@)
	@$(RM) $(TMP_H)

show :
	@echo "MAKEFILE_IN: $(MAKEFILE_IN)"
	@echo "CC: $(CC)"
	@echo "CFLAGS: $(CFLAGS)"
	@echo "EXE: $(EXE)"
	@echo "OBJ: $(OBJ)"
	@echo "SRC: $(SRC1) $(SRC2)"
	@echo "DEP: $(DEP)"
	@echo "#define M2_HAVE_MPI $(HAVE_MPI)"
	@echo "#define M2_HAVE_GLUT $(HAVE_GLUT)"
	@echo "#define M2_GIT_SHA \"$(GIT_SHA)\""
	@echo "#define M2_DEBUG_FILE $(DEBUG_FILE)"
	@echo "#define M2_INSTALL_PATH \"$(PWD)\""

clean :
	$(MAKE) -C lua clean
	$(MAKE) -C zulu clean
	$(MAKE) -C glui clean
	$(RM) $(OBJ) $(DEP) $(EXE)

$(LUA_L) : .FORCE
	$(MAKE) -C lua $(PLAT)

$(ZULU_L) : .FORCE
	$(MAKE) -C zulu lib MAKEFILE_IN=$(MAKEFILE_IN)

$(GLUI_L) : .FORCE
	$(MAKE) -C glui lib MAKEFILE_IN=$(MAKEFILE_IN)

.FORCE :
.PHONY : show clean

-include *.dep
