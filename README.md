# Ram pressure stripping simulations 
## - the co-evolution of interstellar and circumgalactic media
---
:memo: [`./hydro-rps/`](./hydro-rps) is the *codebase* used to simulate our galaxy (with both CGM and ISM) facing *constant ram pressure* from intracluster medium 
using [PLUTO](https://plutocode.ph.unito.it/documentation.html) (version 4.2 patch 2)

> [!IMPORTANT]
> ### Check out my vlog on [the detailed steps to use Catalyst with PLUTO](https://sites.google.com/view/ritalighosh/use-catalyst-with-your-simulations?authuser=0).

### 🏃 Initialization
* Compile Pluto with its dependencies, and create the executable `./pluto`
  ```
  make -j($nproc) && make clean 
  ```
* Use the provided [slurm jobscript](./hydro-rps/slurm-script) to run parallelly on your supercomputer
  ```
  sbatch slurm-script
  ```
 > **_NOTE:_** 
  Setting up the spack environment on clusters might be difficut sometime.
  [Here's](https://sites.google.com/view/ritalighosh/home-spack-environment?authuser=0) a quick documentation I prepared to use on our supercomputer Param Pravega.
  
#### __:bookmark:Important Files:__ ##
> - [init.c](./hydro-rps/init.c) &rarr; problem initialization and on-the-fly analysis
> - [pluto.ini](./hydro-rps/) &rarr; input parameters
> - [definitions.h](./hydro-rps/definitions.h) &rarr; problem / solver settings
> - [userdef_output.c](./hydro-rps/userdef_output.c) &rarr; script for user-defined outputs
> - [*-pipeline.c](./hydro-rps/) &rarr; python files to be used during runtime by Catalyst
> - [generateCatalystAdaptor.py](./hydro-rps/generateCatalystAdaptor.py) &rarr; generates the [`CatalystAdaptor.h`](./hydro-rps/CatalystAdapotor.h) file
