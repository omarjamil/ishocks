This is an internal shocks model for relativistic jets. Please see the following paper for more details: https://arxiv.org/abs/0909.1309


## SECTION 1: Compile, Run and Results

The precompile version should work.

If not, and if on Linux, run build.linux script:
      > build.linux build
      The above will default to release type build i.e. optimized
      Or you can change build type
      > build.linux build debug/release
you can also use build.linux to clean up:
      > build.linux clean

Otherwise or other systems cd to "BUILD" folder and type: 

      > cmake -DCMAKE_BUILD_TYPE=release ..

CMAKE will produce the necessary build file. On linux this will be a Makefile
(requires c/c++ compilers and stdc++ to be installed on your machine): 
    
      > make

This should create an "ishocks" executable in the folder one up from
BUILD. To run please go to the folder above BUILD and type:

      > ./ishocks

That should create four files "resultsFile.dat", "lightcurve.dat" (you
can change their name in parameters.par), "avgSpecInd.dat" and
"paraFile.dat" (you can ignore this one, I used this for param search runs). 

Column 1 in "lightcurve.dat" is time in seconds. 

Column 2 is the radio(or nu_1) flux in Jy

Column 3 is the IR (nu_2) flux in Jy

Column 4 is the spectral index. In this case it is the ratio of nu_2/nu_1 
(higher over lower).


"resultsFile.dat" can be used to plot indvidual shells along the jet.

Column 1 in "resultsFile.dat" is the position along the jet.

Column 2 is the radio(or nu_1) flux in Jy.

Column 3 is the IR (or nu_2) flux in Jy.

Column 4 is the individual shell length

"avgSpecInd.dat"

Column 1 is the average nu_1 flux in Jy.

Column 2 is the average nu_2 flux in Jy.

Column 3 is the average spec index.

Please delete or move the above files before another run otherwise the
new data will be appended to the previous files.

Once it is compiled, you only need to change the "parameters.par" and no need
to recompile. If you do not delete or mv the older resultsFile and lightcurve
files then the new simulation results will be appended to them.

The standard deviation is hard wired in and 
is set to be very low i.e. almost a delta function.

The "doc_parameters.txt" file contains the description of all the parameters.


## SECTION 2: Paramters Combinations: Avoid any pitfalls:


Case 1: Random Distribution of shells:
        In order to have the injected shell parameters (bulk Lorentz
        factor, Mass, shell width (along the jet) etc.), randomly
        sampled from given distributions please make sure the
        following are done in the parameters file:
        
                        use_shell_file = 0
                        use_data_file = 0

        Please refer to "doc_parameters.txt" for more details on the
        other parameters.

Case 2: Shell distribution from a file:
        Can read the shell parameters file. To switch this on have the
        following:
                
                        use_shell_file = 1
                        use_data_file = 0
                        
        The shell file name can be changed:
               
                        shell_file = <filename>
                        
        The shell file needs to contain the following information in
        separate columns:
        column1:  <Shell injection time> 
        column2:  <shell rest mass> 
        column3:  <bulk Lorentz factor>
        column4:  <shell width>
        Note: Please make sure there is NO new line at the end of this file. 
        Note: Please make sure there is !!NO NEW LINE!! at the end of this file.                        
        
Case 3: Sampling the shells at more than "events"*:
        If you wish to sample radiation from the emitting shells at
        times other than just injection or shell merger then switch
        the following on:
        
                        increase_time_resolution = 1
                        
        This still does not mean the time step is completely discrete
        as the shells will still be sampled at injection and
        collisions as well. The increased time step resolution can be
        set by:

                         step_resolution = <time in seconds e.g 1.>
                         
        The total run time of the simulation is set:

                         total_run_duration = <in seconds> 

        The following sets the time of the final shell injection:
         
                         shell_inj_duration = <in seconds>**      
                              
Case 4: Single/Few blob(s) with expansion into vacuum***:
        If you wish only to change the energy densities of the shells
        and no work is being carried out then:

                         in_vacuum_expansion = 1
        
        We are going to assume the shells in this case are injected
        with certain amount of internal energy**** 
        
                         inj_int_energy = 1
                    
        The fraction of rest mass energy released to internal energy:

                         rel_mass_frac = <0 to 1, percent of
                                         relativistic mass>

        It is worth noting that under this regime you wish to give
        shells certain amount of thermal energy then the shells
        injected !!MUST!! have different bulk lorentz factors. This is due
        to the merged shell BLF calculations, which break down if the
        two colliding shells have the same BLF. In the case the
        collision is happening due to longitudinal expansion of the shell.


*an event is either injection of a new shell or collision/merger of
 two shells already in the jet. In order to increase efficiency the
 shells can sampled for radiation at only these events. This means the
 time steps are not uniform. 

**Please note this is important only important when in Case
  1. Otherwise the shellFile determines the time of the final shell
  injection. Also note that number of shells (N)injected when in Case 1
  is determined by this and "avg_ejection_gap".
 
     N ~ shell_inj_duration / avg_ejection_gap

***This assumption means that no work is being done against any
   external gas with the expansion of conical frustums. Under the
   normal circumstances there is an implicit assumption about external
   (the jet) gas pressure leading to a conical shape. Vacuum assumes
   free expansion.

****Under normal circumstances it is assumed the shells are injected
    with all their energy in the form of bulk kinetic energy and the
    internal energy is only generated when two shells collide. The
    internal energy is necessary if shells are going to emit radiation.
    This can be done in other cases as well, simply means the shells
    will be emitting from the moment they are injected.

If a large gap between ejections and high longitudinal expansion rate
then increased sampling maybe necessary to avoid miscalculation of
shell wall positions. Smaller increments leads to higher accuracy for
calculating shell positions. This happens when shells injected with
internal energy and then try to expand in the negative x-direction or
pile up at the source. 
