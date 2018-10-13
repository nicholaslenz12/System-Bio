## The purpose of the directory is to simulate the growth of the two bacteria populations under the influence of some antibiotic concentration A. 

## The pipeline for plotting a stabilizable region is:
* Open StabilizableRegion.m. First you must choose the values of rho and R to simulate. These are defined by rhoSteps, rhoRange, rSteps, rRange.
* If you have no simulation data yet:
    * If you have no datafile to load from, leave fileToLoad as the empty string, "". loadPoint should be 1 in this case.
    * Run StabilizableRegion.m. Periodically, StabilizableRegion.m will output to the commandline how many simulations have been run.
    * If you end the simulation early, write the matrix z to a file using dlmwrite(). Make sure to record the number of simulations that were completed. Now you have a datafile to load from.
* If you have simulation data:
    * Make sure fileToLoad is the name of the datafile to load from.
    * Set loadPoint to the simulation reached on the previous run.
     * If you end the simulation early, write the matrix z to a file using dlmwrite(). Make sure to record the number of simulations that were completed.
