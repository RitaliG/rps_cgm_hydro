## Ram pressure stripping simulations 
### - the co-evolution of interstellar and circumgalactic media
---
*Codebase* to simulate a galaxy facing *constant ram pressure* from intracluster medium 
using ![PLUTO](https://plutocode.ph.unito.it/documentation.html) (version 4.2 patch 2)

### Check out ![the detailed steps to use Catalyst with PLUTO](https://sites.google.com/view/ritalighosh/use-catalyst-with-your-simulations?authuser=0)

### Initialization
* Compile Pluto with its dependencies, and create the executable `./pluto`
  ```
  make -j($nproc) && make clean 
  ```
* Use the provided slurm jobscript to run parallelly on your supercomputer
  ```
  sbatch slurm-script
  ```
## __:bookmark:Important Files:__ ##
> - init.c &rarr; problem initialization and on-the-fly analysis
> - pluto.ini &rarr; input parameters
> - definitions.h &rarr; problem / solver settings
> - userdef_output.c &rarr; script for user-defined outputs
> - *-pipeline.c &rarr; python files to be used during runtime by Catalyst
> - generateCatalystAdaptor.py &rarr; generates the `CatalystAdaptor.h` file
