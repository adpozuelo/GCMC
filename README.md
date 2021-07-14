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

MC is designed to use the massive parallel paradigm for intensive computation. Thus, MC needs a NVIDIA's GPGPU with CUDA arquitecture.

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

<p> Compile: </p>

	cd GCMC/src && make && cd ..

<p> Execute MC application: </p>

	./mc.exe
  ./mc.exe cpu
	./mc.exe gpu

GPU's memory (shared or not shared)
==========

<p> Two versions of GPU's memory developments are supplied: shared and not shared memory. </p>
<p> By default GCMC is configured to use GPU's shared memory and 128 GPU's threads per block. </p>
<p> Next table shows several tests which are been made to compare both options. </p>
<p> Performance study: shared vs not shared GPU memory. </p>
<p>
MonteCarlo NVT Lennard Jones<br>
Shift enabled. gpu[2080Ti] device<br>
500000 particles: 500000O (FCC)<br>
Energy units: LJ<br>
Temperature: 5.0<br>
Volume: 624999.938<br>
ε: 1.0, σ: 1.0, rc: 8.0<br>
Density: 0.8<br>
</p>

<table>
  <tr>
    <th>th</th>
    <th>sha</th>
    <th>no_sha</th>
    <th>s_up</th>
  </tr>
  <tr>
    <td>64</td>
    <td>398.91</td>
    <td>465.45</td>
    <td>1.17</td>
  </tr>
  <tr>
    <td>128</td>
    <td>358.20</td>
    <td>513.78</td>
    <td>1.43</td>
  </tr>
  <tr>
    <td>256</td>
    <td>379.36</td>
    <td>459.50</td>
    <td>1.21</td>
  </tr>
  <tr>
    <td>512</td>
    <td>400.27</td>
    <td>455.84</td>
    <td>1.14</td>
  </tr>
</table>

<p>
th: GPU threads per block. <br>
sha: average time (sec) with GPU shared memory. <br>
no_sha: average time (sec) without GPU shared memory. <br>
sp_up: no_sha / sha. <br>
</p>

<p> To use shared memory: </p>

	cd src && cp gpu_shared.cu gpu.cu && make && cd ..

<p> To don`t use shared memory: </p>

	cd src && cp gpu_no_shared.cu gpu.cu && make && cd ..