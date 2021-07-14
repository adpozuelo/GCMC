Serial (CPU) and parallel (GPGPU/CUDA) general purpose Monte Carlo
code for atomistic simulations.
===========

<p> Developer: Antonio Díaz Pozuelo (adpozuelo@gmail.com) </p>

<p> MC is a Metropolis Monte Carlo which simulates bulk systems composed of soft spherical particles without charges. </p>

<p> There are examples of input data files in 'input_templates' directory. </p>

<p> The results of this software has been validated against LAMMPS simulations. You are free to use lammps input scripts ('lammps' directory) and validations Python files to check it. </p>

- Implemented ensembles:

  * NVT.

- Implemented energy potentials:

  * Lennard Jones.

- Implemented units:

  * Energy: ElectronVolts, Kelvins and Lennard Jones.
  * Distance: Angstroms and Lennard Jones.

Design & Arquitecture
==========

MC is designed to use the massive parallel paradigm for intensive computation.Thus, MC needs a NVIDIA's GPGPU with CUDA arquitecture.

Implementation
==========
MC is implemented with both NVIDIA's CUDA C/C++ and Intel C/C++ programming languages.

Requisites
==========

- Software:

  * NVIDIA CUDA Compiler (nvcc).
  * Intel C++ Compiler (icpc) with MKL (Math Kernel Library) support.

- Hardware:

  * NVIDIA's GPGPU CUDA capable arquitecture with computational capability 6 or higher.

Install
=======

<p> Download MC application: </p>

	git clone https://github.com/adpozuelo/GCMC.git 

<p> Compile</b>: </p>

	cd GCMC/src && make && cd ..

<p> Execute MC application: </p>

	./mc.exe cpu  <- CPU (serial) mode.
	./mc.exe gpu [gpu_id]  <- GPU (parallel) using CUDA device 0 (default).

GPU memory
==========

<p> Two versions of GPU memory developments are supplied: shared and no shared memory. </p>
<p> Next table explains tests made with both options </p>

<p>
MonteCarlo NVT Lennard Jones
Shift enabled. gpu[2080Ti] device
500000 particles: 500000O (FCC)
Energy units: LJ
Temperature: 5.0
Volume: 624999.938
ε: 1, σ: 1, rc: 8.0
Density: 0.8

Performance study: shared vs no shared GPU memory

th    sha  no_sha s_up
64  398.91 465.45 1.17
128 358.20 513.78 1.43
256 379.36 459.50 1.21
512 400.27 455.84 1.14

th: GPU threads per block
sha: average time (sec) with GPU shared memory
no_sha: average time (sec) without GPU shared memory
sp_up: no_sha / sha
</p>
