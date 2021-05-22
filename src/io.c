/* Author: adpozuelo@gmail.com
 * Version: 1.1
 * Date: 05/2021
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io.h"
#include "conf.h"
#include "energy.h"

#define KBEV 8.6173324e-5

void read_input_file(Configuration *cxf)
{
    char line[MAX_LINE_SIZE];
    char buffer[MAX_LINE_SIZE];
    FILE *file = fopen(INPUT_FILENAME, "r");

    fgets(line, MAX_LINE_SIZE, file);
    sscanf(line, "%u %u %hu %hu %s", &cxf->nstep, &cxf->nequil, &cxf->shift,
           &cxf->lammpstrj, buffer);
    cxf->units = (char *)malloc(strlen(buffer) * sizeof(char));
    strcpy(cxf->units, buffer);
    if (strcmp(cxf->units, "K") != 0 && strcmp(cxf->units, "eV") != 0 && strcmp(cxf->units, "LJ") != 0)
    {
        printf("ERROR: '%s' not implemented as energy unit!\n", cxf->units);
        exit(1);
    }

    cxf->rdmax = (precision *)malloc(NDIM * sizeof(precision));
    double acceptance, rdmax0, rdmax1, rdmax2;
    fgets(line, MAX_LINE_SIZE, file);
    sscanf(line, "%lf %lf %lf %lf", &acceptance, &rdmax0, &rdmax1, &rdmax2);
    cxf->acceptance = acceptance;
    cxf->rdmax[0] = (precision)rdmax0;
    cxf->rdmax[1] = (precision)rdmax1;
    cxf->rdmax[2] = (precision)rdmax2;

    double temp, density;
    fgets(line, MAX_LINE_SIZE, file);
    sscanf(line, "%s %hu %lf %lf", buffer, &cxf->nsp, &temp, &density);
    cxf->temp = (precision)temp;
    cxf->density = (precision)density;
    cxf->input_conf = (char *)malloc(strlen(buffer) * sizeof(char));
    strcpy(cxf->input_conf, buffer);
    if (strcmp(cxf->input_conf, "fcc") != 0 && strcmp(cxf->input_conf, "lammps") != 0)
    {
        printf("ERROR: '%s' not implemented as input configuration!\n", cxf->input_conf);
        exit(1);
    }

    if (strcmp(cxf->units, "eV") == 0)
        cxf->kt = KBEV * cxf->temp;
    if (strcmp(cxf->units, "K") == 0 || strcmp(cxf->units, "LJ") == 0)
        cxf->kt = cxf->temp;

    cxf->nitmax = (cxf->nsp * cxf->nsp + cxf->nsp) / 2;
    cxf->nspps = (unsigned short *)malloc(cxf->nsp * sizeof(unsigned short));
    cxf->atoms = (char **)malloc(cxf->nsp * sizeof(char *));
    cxf->natoms = 0;
    for (int i = 0; i < cxf->nsp; ++i)
    {
        fgets(line, MAX_LINE_SIZE, file);
        sscanf(line, "%s %hu", buffer, &cxf->nspps[i]);
        cxf->atoms[i] = (char *)malloc(strlen(buffer) * sizeof(char));
        strcpy(cxf->atoms[i], buffer);
        cxf->natoms += cxf->nspps[i];
    }

    precision **aa = (precision **)malloc(cxf->nsp * sizeof(precision *));
    precision **bb = (precision **)malloc(cxf->nsp * sizeof(precision *));
    cxf->al = (precision *)malloc(cxf->nitmax * sizeof(precision));
    cxf->bl = (precision *)malloc(cxf->nitmax * sizeof(precision));
    cxf->bl2 = (precision *)malloc(cxf->nitmax * sizeof(precision));
    cxf->rc2 = (precision *)malloc(cxf->nitmax * sizeof(precision));
    cxf->rc = (precision **)malloc(cxf->nsp * sizeof(precision *));
    cxf->itp = (unsigned short **)malloc(cxf->nsp * sizeof(unsigned short *));
    for (int i = 0; i < cxf->nsp; ++i)
    {
        aa[i] = (precision *)malloc(cxf->nsp * sizeof(precision));
        bb[i] = (precision *)malloc(cxf->nsp * sizeof(precision));
        cxf->rc[i] = (precision *)malloc(cxf->nsp * sizeof(precision));
        cxf->itp[i] = (unsigned short *)malloc(cxf->nsp * sizeof(unsigned short));
    }

    double aa_d, bb_d, rc;
    int nit = 0;
    for (int i = 0; i < cxf->nsp; ++i)
    {
        for (int j = i; j < cxf->nsp; ++j)
        {
            fgets(line, MAX_LINE_SIZE, file);
            sscanf(line, "%lf %lf %lf", &aa_d, &bb_d, &rc);
            aa[i][j] = (precision)aa_d;
            bb[i][j] = (precision)bb_d;
            cxf->rc[i][j] = (precision)rc;
            cxf->al[nit] = aa[i][j];
            cxf->bl[nit] = bb[i][j];
            cxf->rc[j][i] = cxf->rc[i][j];
            cxf->rc2[nit] = cxf->rc[i][j] * cxf->rc[i][j];
            cxf->itp[i][j] = nit;
            cxf->itp[j][i] = nit;
            ++nit;
        }
    }

    for (int i = 0; i < cxf->nsp; ++i)
    {
        free(aa[i]);
        free(bb[i]);
    }
    free(aa);
    free(bb);
    fclose(file);

    cxf->sigma_o = cxf->bl[0];
    for (int i = 0; i < cxf->nitmax; ++i)
    {
        if (cxf->bl[i] > cxf->sigma_o)
            cxf->sigma_o = cxf->bl[i];
    }

    for (int i = 0; i < cxf->nitmax; ++i)
    {
        cxf->bl[i] /= cxf->sigma_o;
        cxf->bl2[i] = cxf->bl[i] * cxf->bl[i];
        cxf->al[i] /= cxf->kt;
    }

    for (int i = 0; i < cxf->nsp; ++i)
        for (int j = 0; j < cxf->nsp; ++j)
            cxf->rc[i][j] /= cxf->sigma_o;

    nit = 0;
    for (int i = 0; i < cxf->nsp; ++i)
    {
        for (int j = i; j < cxf->nsp; ++j)
        {
            cxf->rc2[nit] = cxf->rc[i][j] * cxf->rc[i][j];
            ++nit;
        }
    }

    cxf->esrrc = (precision *)malloc(cxf->nitmax * sizeof(precision));
    for (int i = 0; i < cxf->nitmax; ++i)
        if (cxf->shift)
            cxf->esrrc[i] = lennard_jones(cxf->rc2[i], i, cxf);
        else
            cxf->esrrc[i] = 0.0;
}

void write_configuration(const Configuration *cxf, const unsigned int timestep)
{
    FILE *fp;

    if (timestep == 0)
    {
        fp = fopen(OUTPUT_XYZ_FILENAME, "w");
        fclose(fp);
    }

    fp = fopen(OUTPUT_XYZ_FILENAME, "a");
    fprintf(fp, "ITEM: TIMESTEP\n%u\n", timestep);
    fprintf(fp, "ITEM: NUMBER OF ATOMS\n%u\n", cxf->natoms);
    fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n");

    for (int i = 0; i < NDIM; ++i)
        fprintf(fp, "%.15le %.15le\n", 0.0, cxf->side[i] * cxf->sigma_o);

    fprintf(fp, "ITEM: ATOMS id type mol x y z\n");
    for (unsigned int i = 0; i < cxf->natoms; ++i)
    {
        fprintf(fp, "%u %hu %hu ", i + 1, cxf->ptype[i] + 1, cxf->molecule[i]);
        for (int j = 0; j < NDIM; ++j)
            fprintf(fp, "%.15lf ", cxf->r[i * NDIM + j] * cxf->side[j] * cxf->sigma_o);
        fputc('\n', fp);
    }

    fclose(fp);
}

void print_header(const Configuration *cxf)
{
    printf("MonteCarlo NVT Lennard Jones\n");
    if (cxf->shift)
        printf("Shift enabled. %s[%d] device\n", cxf->mode, cxf->cuda_device);
    else
        printf("Shift not enabled. %s[%d] device\n", cxf->mode, cxf->cuda_device);

    printf("%d steps [%d equilibrium]\n", cxf->nstep, cxf->nequil);
    printf("%d particles: ", cxf->natoms);
    for (int i = 0; i < cxf->nsp; ++i)
        printf("%d%s ", cxf->nspps[i], cxf->atoms[i]);
    putchar('\n');

    char len_units = 'A';
    char temp_units = 'K';
    if (strcmp(cxf->units, "LJ") == 0)
    {
        len_units = ' ';
        temp_units = ' ';
    }

    printf("Temperature: %.2lf%c. Energy units: %s\nDensity: %lf natoms/%c^3\nVolume: %.3lf %c^3: ",
           cxf->temp, temp_units, cxf->units, cxf->density, len_units,
           cxf->volume * cxf->sigma_o * cxf->sigma_o * cxf->sigma_o, len_units);
    for (int i = 0; i < NDIM; ++i)
        printf("%.4lf%c ", cxf->side[i] * cxf->sigma_o, len_units);
    putchar('\n');
    printf("--------------------------------------------------\n");
    printf("%-10s  %-6.2s %-22.8s\n", "Step", "% Acc", "PotEng");
    printf("--------------------------------------------------\n");

    FILE *fp = fopen(LOG_FILENAME, "w");
    fputs("MonteCarlo NVT Lennard Jones\n", fp);
    if (cxf->shift)
        fprintf(fp, "Shift enabled. %s[%d] device\n", cxf->mode, cxf->cuda_device);
    else
        fprintf(fp, "Shift not enabled. %s[%d] device\n", cxf->mode, cxf->cuda_device);

    fprintf(fp, "%d steps [%d equilibrium]\n", cxf->nstep, cxf->nequil);
    fprintf(fp, "%d particles: ", cxf->natoms);
    for (int i = 0; i < cxf->nsp; ++i)
        fprintf(fp, "%d%s ", cxf->nspps[i], cxf->atoms[i]);
    fputc('\n', fp);

    fprintf(fp, "Temperature: %.2lf%c. Energy units: %s\nDensity: %lf natoms/%c^3\nVolume: %.3lf %c^3: ",
           cxf->temp, temp_units, cxf->units, cxf->density, len_units,
           cxf->volume * cxf->sigma_o * cxf->sigma_o * cxf->sigma_o, len_units);
    for (int i = 0; i < NDIM; ++i)
        fprintf(fp, "%.4lf%c ", cxf->side[i] * cxf->sigma_o, len_units);
    fputc('\n', fp);
    fputs("--------------------------------------------------\n", fp);
    fprintf(fp, "%-10s  %-6.2s %-22.8s\n", "Step", "% Acc", "PotEng");
    fputs("--------------------------------------------------\n", fp);
    fclose(fp);
}

void print_step(const Configuration *cxf, const unsigned int step)
{
    precision naccepter = (cxf->naccept / (precision)cxf->ntrial) * 100.0;
    printf("%-10u  %-6.2lf %-22.8lf\n", step, naccepter, cxf->esr * cxf->kt);

    FILE *fp = fopen(LOG_FILENAME, "a");
    fprintf(fp, "%-10u  %-6.2lf %-22.8lf\n", step, naccepter, cxf->esr * cxf->kt);
    fclose(fp);
}