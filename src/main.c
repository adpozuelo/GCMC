/*
 * Serial (CPU) and parallel (GPGPU/CUDA) general purpose Monte Carlo
 * code for atomistic simulations.
 *
 * At present simulates bulk systems composed of soft spherical particles
 * without charges.
 * 
 * Implemented ensembles: NVT.
 * Implemented energy potentials: Lennard Jones.
 * 
 * Program units:
 *        Energy: "eV" electronVolts.
 *                 "K" Kelvins.
 *        Distance: Angstrom.
 *
 * Input data files (see the files for parameter specifications):
 *          in.mc : contains the configuration of the system to be simulated.
 *          init.lammpstrj : initial configuration (only for lammps option in in.mc file)
 *
 * Output files:
 *           log.mc : simulation output
 *           mc.lammpstrj : trajectory file
 *
 * Single or double precision can be set in conf.h header.
 * 
 * Author: adpozuelo@gmail.com
 * Version: 1.1
 * Date: 05/2021
 */

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mkl_vsl.h"
#include "conf.h"
#include "io.h"
#include "xyz.h"
#include "energy.h"
#include "moves.h"
#include "gpu.h"

int main(int argc, char *argv[])
{
    if (argc < 2 || argc > 3)
    {
        puts("Usage: mc cpu/gpu [gpuId]");
        exit(1);
    }

    Configuration cxf;

    if (strcmp(argv[1], "cpu") != 0 && strcmp(argv[1], "gpu") != 0)
    {
        puts("Usage: mc cpu/gpu [gpuId] input_file");
        exit(1);
    }
    else
    {
        cxf.mode = (char *)malloc(strlen(argv[1]) * sizeof(char));
        strcpy(cxf.mode, argv[1]);
    }

    cxf.cuda_device = cxf.ntrial = cxf.naccept = 0;
    if (argv[2] != NULL)
        cxf.cuda_device = atoi(argv[2]);

    // vslNewStream(&run.streamRNG, VSL_BRNG_MT19937, time(NULL)); // Random mode
    vslNewStream(&cxf.streamRNG, VSL_BRNG_MT19937, 1); // Determinist mode

    read_input_file(&cxf);

    if (strcmp(cxf.input_conf, "fcc") == 0)
        init_fcc(&cxf);
    if (strcmp(cxf.input_conf, "lammps") == 0)
        read_lammpstrj(&cxf);
    
    if (strcmp(cxf.mode, "gpu") == 0)
        gpu(&cxf, 0);

    print_header(&cxf);

    if (strcmp(cxf.mode, "cpu") == 0)
        cxf.esr = energy_cpu(&cxf);
    if (strcmp(cxf.mode, "gpu") == 0)
        gpu(&cxf, 1);

    print_step(&cxf, 0);
    if (cxf.lammpstrj == 1 || cxf.lammpstrj == 2)
        write_configuration(&cxf, 0);

    for (unsigned int step = 1; step <= cxf.nstep; ++step)
    {
        if (strcmp(cxf.mode, "cpu") == 0)
            move_atoms_cpu(&cxf);
        if (strcmp(cxf.mode, "gpu") == 0)
            gpu(&cxf, 2);

        if (step <= cxf.nequil)
            for (int i = 0; i < NDIM; ++i)
                cxf.rdmax[i] *= (cxf.naccept / (precision)cxf.ntrial) / cxf.acceptance;

        print_step(&cxf, step);
        
        if (step > cxf.nequil && cxf.lammpstrj == 2)
            write_configuration(&cxf, step);
    }

    if (strcmp(cxf.mode, "gpu") == 0)
        gpu(&cxf, 3);
    for (int i = 0; i < cxf.nsp; ++i)
    {
        free(cxf.atoms[i]);
        free(cxf.itp[i]);
        free(cxf.rc[i]);
    }
    vslDeleteStream(&cxf.streamRNG);
    free(cxf.units);
    free(cxf.mode);
    free(cxf.al);
    free(cxf.bl);
    free(cxf.bl2);
    free(cxf.rc);
    free(cxf.rc2);
    free(cxf.itp);
    free(cxf.atoms);
    free(cxf.rdmax);
    free(cxf.nspps);
    free(cxf.r);
    free(cxf.side);
    free(cxf.esrrc);
    free(cxf.input_conf);
    free(cxf.ptype);

    putchar('\n');
    return 0;
}
