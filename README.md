Serial (CPU) and parallel (GPGPU/CUDA) general purpose Monte Carlo
code for atomistic simulations. (1.0 stable version)
===========

<p> Developer: Antonio Díaz Pozuelo (adpozuelo@gmail.com) </p>

<p> MC is a Metropolis Monte Carlo which simulates bulk systems composed of soft spherical particles without charges. </p>

<p> Implemented ensembles: NVT. </p>
<p> Implemented energy potentials: Lennard Jones. </p>

<p> Implemented units:
	- Energy: 
		* ElectronVolts.
 		* Kelvins.
	- Distance:
		* Angstroms
</p>

<p> There are examples of input data files in 'input_templates' directory. </p>

<p> The results of this software has been validated against LAMMPS simulations. You are free to use lammps input scripts ('lammps' directory) and validations Python files to check it. </p>

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

	git clone -- 

<p> Compile </b>: </p>

	cd GCMC/src && make && cd ..

<p> Execute MC application: </p>

	./mc.exe cpu  <- CPU (serial) mode.
	./mc.exe gpu [gpu_id]  <- GPU (parallel) using CUDA device 0 (default).
