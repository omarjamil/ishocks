##########################################################
# Makefile
#
SOURCES.cc = numerical.cc parameters.cc paraccess.cc shells.cc radiation.cc container.cc mergers.cc evolution.cc simulation.cc run.cc
INCLUDES = 
CPP = g++
CFLAGS = 
CPPFLAGS = -O3 -W
LIBS = 

# the executable file
EFILE= ../ishocks
# Object file
OBJECTS= $(SOURCES.cc:.cc=.o)

all: $(EFILE)


# linking to create the executable file
$(EFILE): $(OBJECTS)
	@echo "linking...." 
	$(CPP) $(CPPFLAGS) -o $(EFILE) $(OBJECTS) $(LIBS)


# this declares dependence for object files
# so only the necessary ones
# if the dependence files changed, the .o file
# will be recompiled

numerical.o: numerical.hh numerical.cc
	$(CPP) $(CPPFLAGS) -c numerical.cc

parameters.o: parameters.hh parameters.cc 
	$(CPP) $(CPPFLAGS) -c parameters.cc 

paraccess.o: paraccess.cc paraccess.hh
	$(CPP) $(CPPFLAGS) -c paraccess.cc

shells.o: shells.hh shells.cc paraccess.hh
	$(CPP) $(CPPFLAGS) -c shells.cc

radiation.o: radiation.hh radiation.cc paraccess.hh shells.hh shells.cc container.hh container.cc 
	$(CPP) $(CPPFLAGS) -c radiation.cc

container.o: container.hh container.cc shells.hh shells.cc paraccess.hh numerical.hh
	$(CPP) $(CPPFLAGS) -c container.cc

mergers.o: mergers.hh mergers.cc shells.hh container.hh container.cc funcobj.hh paraccess.hh
	$(CPP) $(CPPFLAGS) -c mergers.cc

evolution.o: evolution.hh evolution.cc mergers.cc funcobj.hh simulation.cc paraccess.hh
	$(CPP) $(CPPFLAGS) -c evolution.cc

simulation.o: simulation.hh simulation.cc funcobj.hh shells.hh paraccess.hh
	$(CPP) $(CPPFLAGS) -c simulation.cc

run.o: run.cc simulation.hh simulation.cc 
	$(CPP) $(CPPFLAGS) -c run.cc

##########################################################
# run "make debug"
# "./debug" will produce "gmon.out"
# "gprof debug" can be used to analyze the output
debug:
	g++ -pg $(OBJECTS) $(LIBS) -o ../$@


# remove the object and executable files
clean:
	rm -f $(EFILE) $(OBJECTS)

spotless:
	rm -f $(EFILE) $(OBJECTS) *~

#
##########################################################