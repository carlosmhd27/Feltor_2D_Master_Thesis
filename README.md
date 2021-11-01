# FELTOR 2D

Repository for the Master's thesis in [_Transition from drift-wave to interchange turbulence in magnetically confined plasmas_](https://www.overleaf.com/read/qqdnjjzyrwnv) by Carlos Rodríguez.

___
## Abstract

Throughout this work, we present a novel technique to characterize drift-wave and interchange turbulence in the limits of the last closed flux surface. This method is based on the study of the phase relation between electrostatic potential and density as well as between radial velocity and density using the cross power spectral density. 2D simulations, written using [FELTOR](https://feltor-dev.github.io/), of an interchange and the modified Hasegawa-Wakatani models are studied. Using these simulations we are able to measure and provide the differences between drift-wave and interchange instability. When measuring the phase shift between potential and density, we get a 0 rad difference for drift-wave turbulence and <img src="https://latex.codecogs.com/svg.image?\frac{\pi}{2}" title="\frac{\pi}{2}" /> rad for interchange. When comparing the phase difference between radial transport and density, the results are the opposite, 0 rad for interchange and <img src="https://latex.codecogs.com/svg.image?\frac{\pi}{2}" title="\frac{\pi}{2}" /> rad for drift-wave. Afterwards, a complete model is proposed where plasma dynamic is studied in the regions around the last closed flux surface, the edge and the scrape-off layer. From these simulations we measure a mixture of interchange and drift-wave turbulence in the limits of the edge region, where the phase relation between radial velocity and density becomes <img src="https://latex.codecogs.com/svg.image?\frac{\pi}{4}" title="\frac{\pi}{4}" /> rad. This allow us to measure the contribution of both instabilities when they co-exist. Furthermore, we find a change of regime when crossing the last closed flux surface, the phase relation drops to 0 rad at the last closed flux surface, which can be used to measure the position of it. Finally, this technique is applied to state-of-the-art experimental data from the EAST tokamak reactor, where different phase relations are measured at the edge and at the scrape-off layer. The phase relation between electron temperature and floating potential change from _π_ rad in the edge to <img src="https://latex.codecogs.com/svg.image?-\frac{3\pi}{4}" title="-\frac{3\pi}{4}" /> rad in the scrape-off layer. This result could be used to locate the last closed flux surface, reducing the cost of this process.
___

During this work 4 models where implemented to the study of plasma dynamics and turbulence inside a magnetically confinement nuclear fusion reactor, in the boundaries of the Last Closed Flux Surface (LCFS). This models separate the study of Drift-Wave turbulence, in the inner region of the reactor, by using the Hasegawa-Wakatani model [[1]](#HW) [[2]](#HW1) --in the ordinary and modified versions--; the study of Interchange turbulence in the Scrape-Off Layer (SOL), using an Interchange model  [[3]](#IC); and the study of both turbulence interacting together, by using a complete model [[4]](#HWIC). The partial differential equations for these models can be seen in the documentation, `docuemntation.pdf`. The equations are solved using Discontinuous Galerkin methods with the [FELTOR](https://feltor-dev.github.io/) scientific library in C++. To install and compile the necessary dependencies, one should follow the [FELTOR documentation in Github](https://github.com/feltor-dev/feltor). In general, you need to clone the dependencies of [FELTOR](https://github.com/feltor-dev/feltor), [thrust](https://www.github.com/thrust/thrust), [cusp](https://www.github.com/cusplibrary/cusplibrary) and [vlc](https://www.github.com/vectorclass/version1)

```
git clone https://www.github.com/feltor-dev/feltor
git clone https://www.github.com/thrust/thrust
git clone https://www.github.com/cusplibrary/cusplibrary
git clone https://www.github.com/vectorclass/version1 vcl
```

Thuse of [FELTOR](https://feltor-dev.github.io/) allows an easy optimization of the code using parallelization techniques, vectorization, multithreading and MPI, as well as GPU optimization using [CUDA](https://developer.nvidia.com/cuda-downloads). To compile any of this alternatives in the different directories, it is simply done by using the make file.

```
make convection_hpc device=omp #(for an OpenMP version)
#or
make convection_hpc device=gpu #(if you have a gpu and nvcc )
#or
make convection_mpi device=<...> #(if you want to run an MPI version)

```

The MPI version is available for the `FELTOR-MPI` and `FELTOR-Complete` models.

All the directories present 2 compilation forms, `convection` and `convection_hpc` (or `convection_mpi`). The first one allows to visualize the evolution of the initial condition, using [GLFW3](https://www.glfw.org/) while the second one saves the output in a [netCDF4](https://www.unidata.ucar.edu/software/netcdf/) file. In both cases the initial conditions are given through a [JSON](https://www.json.org/json-en.html) file. Hence to run the binary programs, one must run:

```
./convection input.json #(for visualizing)
#or
./convection_hpc input.json output.nc #(for saving)
```

To run the mpi version, one has to select the number of processors and the distribution of them in a square grid:

```
mpirun -n 6 ./convection_mpi input.json output.nc 2 3 #(for example)
```

The intial conditions might vary from model to model, read the documentation `documentation.pdf` for more information, every folder has already a working initial condition file.

The Hasegawa - Wakatani and Interchange models are run separately with the models in the `FELTOR-model` and `FELTOR-MPI` folders. The Complete model is run in the `FELTOR-Complete` model.

To analyse and animate the output, the `analysis` Python class is created, in `2D_FELTOR_analysis`; and the `GIF_main.py` script in `Python_GIF`. The analysis is done using either probes or the whole grid. It is important to notice that the grid is reduced every certain number of steps to be saved and reduce the weight of the output folder, hence the probes are used to study the evolution of the system in a certain position of the grid in every time steps. For further information and view of the results, one can read the Master's thesis [_Transition from drift-wave to interchange turbulence in magnetically confined plasmas_](https://www.overleaf.com/read/qqdnjjzyrwnv).


Finally, the `Miscellanious` folder has some random python and slurm scripts used to run the program in the [marconi](https://www.hpc.cineca.it/hardware/marconi) and [marconi100](https://www.hpc.cineca.it/hardware/marconi100) high performance computing CPU and GPU clusters from [CINECA](https://www.cineca.it/en/), as this work has been part of the [EUROfusion](https://www.euro-fusion.org/) organization to contribute in the reasearch of Nuclear Fusion.



## References
<a id="HW">[1]</a>
S. B. Korsholm, _Coherent structures and transport in drift wave plasma turbulence_, Ph.D. dissertation, Danmarks Tekniske Universitet (DTU), 2011.

<a id="HW2">[2]</a>
J.  Anderson  and  B.  Hnat, _Statistical  analysis  of  hasegawa-wakatani  turbulence_,Physics of Plasmas, vol. 24, p. 062 301, 2017. DOI: [10.1063/1.4984985](https://doi.org/10.1063/1.4984985).

<a id="IC">[3]</a>
R.  Kube,  O.  E.  García,  and  M.  Wiesenberger, _Amplitude  and  size  scaling  for  interchangemotions of plasma filaments_, Physics of Plasmas, vol. 23, no. 12, p. 122 302, 2016. DOI: [10.1063/1.4971220]((https://doi.org/10.1063/1.4971220)).

<a id="HWIC">[4]</a>
 G. Decristoforo, A. Theodorsen, J. Omotani, T. Nicholas, and O. E. Garcia, _Numerical tur-bulence simulations of intermittent fluctuations in the scrape-off layer of magnetized plasmas_, 2021. arXiv: [2102.04723 [physics.plasm-ph]](https://arxiv.org/abs/2102.04723).
