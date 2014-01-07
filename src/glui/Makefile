
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


# glui executables
# --------------------------------------------------
EXE_SRC = $(wildcard example/*.cpp)
EXE_OBJ = $(EXE_SRC:.cpp=.o)
EXE     = $(EXE_SRC:.cpp=.exe)


# object code required for executables
# --------------------------------------------------
SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)
DEP = $(SRC:.cpp=.dep)
LIB = libglui.a

# build rules
# --------------------------------------------------
default : $(EXE)
lib : $(LIB)

$(LIB) : $(OBJ)
	$(AR) $@ $?
	$(RANLIB) $@

example/%.o : example/%.cpp $(MAKEFILE_IN)
	@$(CXX) -MM $< > $(<:.cpp=.dep) $(GLUT_I) -I$(PWD)
	$(CXX) $(CFLAGS) $< -c -o $@ $(GLUT_I) -I$(PWD)

%.o : %.cpp $(MAKEFILE_IN)
	@$(CXX) -MM $< > $(<:.cpp=.dep) $(GLUT_I)
	$(CXX) $(CFLAGS) $< -c $(GLUT_I)

%.exe : %.o $(OBJ)
	$(CXX) $(CFLAGS) $^ $(CLIBS) $(GLUT_L) -o $@

show :
	@echo "MAKEFILE_IN: $(MAKEFILE_IN)"
	@echo "CC: $(CC)"
	@echo "CFLAGS: $(CFLAGS)"
	@echo "EXE: $(EXE)"
	@echo "OBJ: $(OBJ) $(EXE_OBJ)"
	@echo "SRC: $(SRC) $(EXE_SRC)"
	@echo "DEP: $(DEP)"

clean :
	$(RM) $(OBJ) $(EXE_OBJ) $(DEP) $(EXE) $(LIB)

.FORCE :

-include *.dep
