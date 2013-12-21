
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
ifeq ($(HAVE_GLUT), 1)

endif


# object code required for executables
# --------------------------------------------------
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)
DEP = $(SRC:.c=.dep)


# build rules
# --------------------------------------------------
default : $(EXE)

%.o : %.c m2-cfg.h $(MAKEFILE_IN)
	@$(CC) -MM $< > $(<:.c=.dep)
	$(CC) $(CFLAGS) $< -c

m2 : $(OBJ)
	$(CC) $(CFLAGS) $^ $(CLIBS) -o $@

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
	@echo "SRC: $(SRC)"
	@echo "DEP: $(DEP)"

clean :
	$(RM) $(OBJ) $(DEP) $(EXE)

.FORCE :

-include *.dep
