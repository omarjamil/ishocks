##########################################################
# Makefile
#
SOURCES.cc = statAnalysis.cc run.cc
INCLUDES = 
CC = g++
CFLAGS = -O3
CPPFLAGS = 
LIBS = 

# the executable file
EFILE= ../periodogram
# Object file
OBJECTS= $(SOURCES.cc:.cc=.o)

all: $(EFILE)


# linking to create the executable file
$(EFILE): $(OBJECTS)
	@echo "linking...." 
	$(CC) $(CFLAGS) -o $(EFILE) $(OBJECTS) $(LIBS)


# this declares dependence for object files
# so only the necessary ones
# if the dependence files changed, the .o file
# will be recompiled

statAnalysis.o: statAnalysis.hh statAnalysis.cc
	$(CC) $(CFLAGS) -c statAnalysis.cc

run.o: run.cc statAnalysis.hh statAnalysis.cc 
	$(CC) $(CFLAGS) -c run.cc



# run "make debug"
# "./debug" will produce "gmon.out"
# "gprof debug" can be used to analyze the output
debug:
	g++ -pg $(OBJECTS) $(LIBS) -o $@


# remove the object and executable files
clean:
	rm -f $(EFILE) $(OBJECTS)

spotless:
	rm -f $(EFILE) $(OBJECTS) *~

#
##########################################################