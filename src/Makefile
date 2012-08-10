CC = gcc
# Compiler flags
CFLAGS = # -c -O3 -Wall -g
# C++ compiler flags
CPPFLAGS = -c -O3 -std=c++0x -Wall -g
# Header files
INCLUDE = -I../include/
# Linker flags
LDFLAGS = -g -Wall
# Linker libraries
LDLIBS = -lstdc++ -lm
# Sources
SRCS = $(wildcard *.cc)
# Headers
HDRS = $(wildcard ../include/*.hh)
# the executable file
EFILE= ../ishocks

# Object file
OBJECTS= $(SRCS:.cc=.o)

# linking to create the executable file
# $@ is the name of the target
$(EFILE): $(OBJECTS)
	@echo "linking...." 
	$(CC) $(OBJECTS) $(LDLIBS) $(LDFLAGS) -o $@ 

# $< list first dependency, therefore the name of the file from include depends
# percent sign does text substituion so can be used to write generic makefiles
%.o: %.cc
	$(CC) $(CFLAGS) $(CPPFLAGS) $(INCLUDE) -c $<

clean:
	rm -f $(EFILE) $(OBJECTS)

spotless:
	rm -f $(EFILE) $(OBJECTS) *~

# this generates the dependencies list stored in depend file
# depend is regenerated everytime SRC or HDRS change
depend: $(SRCS) $(HDRS)
	$(CC) -MM $(CPFLAGS) $(INCLUDE) $(SRCS) > depend


# run "make debug"
# "./debug" will produce "gmon.out"
# "gprof debug" can be used to analyze the output
debug:
	$(CC) -pg $(OBJECTS) $(LDLIBS) $(LDFLAGS) -o ../$@

#depend: 
#	$(CC) -MM $(CPFLAGS) $(INCLUDE) $(SRCS) > depend



# include the dependencies file
include depend
