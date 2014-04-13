# --------------------------------------------------
# if there is no Makefile.in then use the template
# --------------------------------------------------
ifneq ($(strip $(MAKEFILE_IN)),)
# use value of MAKEFILE_IN if provided on the command line
else ifeq ($(shell test -e Makefile.in && echo 1), 1)
MAKEFILE_IN = Makefile.in
else
MAKEFILE_IN = Makefile.in.template
endif
include $(MAKEFILE_IN)

B2C = bin2c.lua
INC += -I$(LUA_INC)

ifeq ($(strip $(HAVE_HDF5)), 1)
INC += -I$(HDF5_HOME)/include
endif

SRC = $(wildcard *.c)
DEP = $(SRC:.c=.dep)
ZULU_LIB = libzulu.a
ZULU_OBJ = \
	30log.o \
	ansicolors.o \
	argparse.o \
	array.o \
	buffer.o \
	class.o \
	serpent.o \
	struct.o \
	H5.o \
	hdf5.o \
	mpi.o \
	json.o \
	lfs.o \
	zulu.o

default : $(ZULU_LIB)

lib : $(ZULU_LIB)

show :
	@echo "MAKEFILE_IN=$(MAKEFILE_IN)"
	@echo "HAVE_MPI=$(HAVE_MPI)"
	@echo "HAVE_HDF5=$(HAVE_HDF5)"
	@echo "ZULU_LIB=$(ZULU_LIB)"
	@echo "INC=$(INC)"
	@echo "DEP=$(DEP)"

$(ZULU_LIB) : $(ZULU_OBJ)
	$(AR) $@ $?
	$(RANLIB) $@

%.c : %.lua
	$(LUA_EXE) $(B2C) +$< > $@

%.o : %.c $(MAKEFILE_IN)
	@$(CC) -MM $< $(INC) > $(<:.c=.dep)
	$(CC) $(CFLAGS) -o $@ $< $(INC) -c

clean :
	$(RM) $(ZULU_LIB) $(ZULU_OBJ) $(wildcard *.dep)
