
# if there is no Makefile.in then use the template
# --------------------------------------------------
ifeq ($(shell test -e Makefile.in && echo 1), 1)
MAKEFILE_IN = Makefile.in
else
MAKEFILE_IN =  Makefile.in.template
endif
include $(MAKEFILE_IN)


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
GLUI_L = glui/libglui.a


# build rules
# --------------------------------------------------
default : $(EXE)

vis.o : vis.cpp m2-cfg.h $(MAKEFILE_IN)
	@$(CC) -MM $< > $(<:.cpp=.dep) $(GLUT_I)
	$(CC) $(CFLAGS) $< -c $(GLUT_I)

%.o : %.c m2-cfg.h $(MAKEFILE_IN)
	@$(CC) -MM $< > $(<:.c=.dep)
	$(CC) $(CFLAGS) $< -c

m2 : $(OBJ) $(GLUI_L)
	$(MAKE) -C glui lib
	$(CXX) $(CFLAGS) $(OBJ) $(CLIBS) $(GLUT_L) $(GLUI_L) -o $@

m2-cfg.h : .FORCE
	@echo "/* m2 config header file */" > $(TMP_H)
	@echo "#define M2_HAVE_MPI $(HAVE_MPI)" >> $(TMP_H)
	@echo "#define M2_GIT_SHA \"$(GIT_SHA)\"" >> $(TMP_H)
	@echo "#define M2_DEBUG_FILE $(DEBUG_FILE)" >> $(TMP_H)
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

clean :
	$(MAKE) -C glui clean
	$(RM) $(OBJ) $(DEP) $(EXE)

$(GLUI_L) : .FORCE
	$(MAKE) -C glui lib

.FORCE :
.PHONY : show clean

-include *.dep
