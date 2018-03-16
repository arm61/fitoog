//
// Created by arm61 on 14/01/18.
//

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define PI 3.141592653
#define PHI 1.618033989

struct PosAng
{
    float position[3];
    float angle[3];
};

struct CharPair
{
    char label[50];
    char keyword[50];
};

struct Atom
{
    char label[3];
    float x;
    float y;
    float z;
    float sl[10];
};

struct Job
{
    char name[50];
    int molecule_types_number;
    int restart;
    float cell[3];
    int scattering_data_number;
    int golden_vectors;
    int population_size;
    int steps_number;
    int population_per_core;
    int max_data_length;
    int max_mol_length;
    int max_mol_num;
    float phi[2];
    int print_freq;
};

struct Data
{
    float q;
    float i;
    float di;
    float dq;
};

void check_file_exists(FILE *file, char *fileName)
{
    if (file == NULL)
    {
        printf("The %s file could not be found, as a result fitoog has quit.\n", fileName);
        exit(EXIT_FAILURE);
    }
}

char *find_word_in_line(char line[100], char word[50])
{
    char *found;
    found = strstr(line, word);
    return found;
}

char *first_char_after_space(char *input)
{
    char *starting = input;
    while (*starting != ';')
    {
        starting++;
    }
    starting++;
    int i = 0;
    while (starting[i] != '\0')
    {
        i++;
    }
    starting[i-1] = '\0';
    return starting;
}

void read_input_file(struct Job *job)
{
    FILE *input_file;
    input_file = fopen("fitoog.inp", "r");
    check_file_exists(input_file, "fitoog.inp");
    char line[100];
    char *found;
    char temp[50];
    const char wipe[50] = "                                                  ";
    job->golden_vectors = 30;
    job->restart = 0;
    job->steps_number = 10;
    while(fgets(line, sizeof(line), input_file))
    {
        found = find_word_in_line(line, "name");
        if(found != NULL)
        {
            strncpy(job->name, wipe, 50);
            strncpy(job->name, first_char_after_space(line), 50);
            strcat(job->name, "\0");
        }
        found = find_word_in_line(line, "num_molecules");
        if(found != NULL)
        {
            strncpy(temp, wipe, 50);
            strncpy(temp, first_char_after_space(line), 50);
            strcat(temp, "\0");
            job->molecule_types_number = atoi(temp);
        }
        found = find_word_in_line(line, "phi0");
        if(found != NULL)
        {
            strncpy(temp, wipe, 50);
            strncpy(temp, first_char_after_space(line), 50);
            strcat(temp, "\0");
            job->phi[0] = atof(temp);
        }
        found = find_word_in_line(line, "phi1");
        if(found != NULL)
        {
            strncpy(temp, wipe, 50);
            strncpy(temp, first_char_after_space(line), 50);
            strcat(temp, "\0");
            job->phi[1] = atof(temp);
        }
        found = find_word_in_line(line, "cell0");
        if(found != NULL)
        {
            strncpy(temp, wipe, 50);
            strncpy(temp, first_char_after_space(line), 50);
            strcat(temp, "\0");
            job->cell[0] = atof(temp);
        }
        found = find_word_in_line(line, "cell1");
        if(found != NULL)
        {
            strncpy(temp, wipe, 50);
            strncpy(temp, first_char_after_space(line), 50);
            strcat(temp, "\0");
            job->cell[1] = atof(temp);
        }
        found = find_word_in_line(line, "cell2");
        if(found != NULL)
        {
            strncpy(temp, wipe, 50);
            strncpy(temp, first_char_after_space(line), 50);
            strcat(temp, "\0");
            job->cell[2] = atof(temp);
        }
        found = find_word_in_line(line, "num_data");
        if(found != NULL)
        {
            strncpy(temp, wipe, 50);
            strncpy(temp, first_char_after_space(line), 50);
            strcat(temp, "\0");
            job->scattering_data_number = atoi(temp);
        }
        found = find_word_in_line(line, "pop_size");
        if(found != NULL)
        {
            strncpy(temp, wipe, 50);
            strncpy(temp, first_char_after_space(line), 50);
            strcat(temp, "\0");
            job->population_size = atoi(temp);
        }
        found = find_word_in_line(line, "num_steps");
        if(found != NULL)
        {
            strncpy(temp, wipe, 50);
            strncpy(temp, first_char_after_space(line), 50);
            strcat(temp, "\0");
            job->steps_number = atoi(temp);
        }
        found = find_word_in_line(line, "restart;yes");
        if(found != NULL)
        {
            job->restart = 1;
        }
    }
    fclose(input_file);
}

void get_iterator_labels(int num, struct CharPair pair[num], char id[50])
{
    int i;
    for(i = 0; i < num; i++)
    {
        strncpy(pair[i].label, id, 50);
        char str[50];
        sprintf(str, "%d", i);
        strcat(pair[i].label, str);
        strcat(pair[i].label, "\0");
    }
    FILE *input_file;
    input_file = fopen("fitoog.inp", "r");
    check_file_exists(input_file, "fitoog.inp");
    char line[100];
    while(fgets(line, sizeof(line), input_file))
    {
        int i;
        for (i = 0; i < num; i++)
        {
            char *found = find_word_in_line(line, pair[i].label);
            if (found != NULL)
            {
                const char wipe[50] = "                                                  ";
                strncpy(pair[i].keyword, wipe, 50);
                strncpy(pair[i].keyword, first_char_after_space(line), 50);
                strcat(pair[i].keyword, "\0");
            }
        }
    }
    fclose(input_file);
}


void get_files_lengths(int num_files, int lengths[num_files], struct CharPair pair[num_files])
{
    int i;
    for (i = 0; i < num_files; i++)
    {
        FILE *input_file;
        input_file = fopen(pair[i].keyword, "r");
        check_file_exists(input_file, pair[i].keyword);
        char line[100];
        int j = 0;
        while (fgets(line, sizeof(line), input_file))
        {
            j++;
        }
        lengths[i] = j;
        fclose(input_file);
    }
}

int find_max_int(int length, int array[length])
{
    int i;
    int max = 0;
    for (i = 0; i < length; i++)
    {
        if (array[i] > max)
        {
            max = array[i];
        }
    }
    return max;
}

int find_max_charpair(int length, struct CharPair array[length])
{
    int i;
    int max = 0;
    for (i = 0; i < length; i++)
    {
        float num = atof(array[i].keyword);
        if (num > max)
        {
            max = num;
        }
    }
    return max;
}

void get_data_points(struct Job job, struct Data exp_data[job.scattering_data_number][job.max_data_length],
                     struct CharPair data_types[job.scattering_data_number])
{
    int i;
    for(i = 0; i < job.scattering_data_number; i++)
    {
        FILE *data_file;
        data_file = fopen(data_types[i].keyword, "r");
        check_file_exists(data_file, data_types[i].keyword);
        char line[512];
        int k = 0;
        while(fgets(line, sizeof(line), data_file))
        {
            int j = 0;
            int maxj = 0;
            char *str;
            str = strtok(line, " ");
            while (str != NULL)
            {
                if (j == 0)
                {
                    exp_data[i][k].q = atof(str);
                }
                if (j == 1)
                {
                    exp_data[i][k].i = atof(str);
                }
                if (j == 2)
                {
                    exp_data[i][k].di = atof(str);
                }
                if (j == 3)
                {
                    exp_data[i][k].dq = atof(str);
                }
                maxj = j;
                j++;
                str = strtok(NULL, ",");
            }
            if (maxj < 2)
            {
                exp_data[i][k].di = exp_data[i][k].i * 0.05;
                exp_data[i][k].dq = exp_data[i][k].q * 0.05;
            }
            if (maxj < 3)
            {
                exp_data[i][k].dq = exp_data[i][k].q * 0.05;
            }
            k++;
        }
        fclose(data_file);
    }
}

void get_atom_positions(struct Job job, struct Atom molecules[job.molecule_types_number][job.max_mol_length],
                        struct CharPair mol_types[job.molecule_types_number])
{
    int i;
    for(i = 0; i < job.molecule_types_number; i++)
    {
        FILE *mol_file;
        mol_file = fopen(mol_types[i].keyword, "r");
        check_file_exists(mol_file, mol_types[i].keyword);
        char line[512];
        int k = 0;
        while(fgets(line, sizeof(line), mol_file))
        {
            int j = 0;
            char *str;
            str = strtok(line, " ");
            while (str != NULL)
            {
                if (j == 0)
                {
                    strncpy(molecules[i][k].label, str, 3);
                }
                if (j == 1)
                {
                    molecules[i][k].x = atof(str);
                }
                if (j == 2)
                {
                    molecules[i][k].y = atof(str);
                }
                if (j == 3)
                {
                    molecules[i][k].z = atof(str);
                }
                int r;
                for (r = 0; r < job.scattering_data_number; r++)
                {
                    if (j == r + 4)
                    {
                        molecules[i][k].sl[r] = atof(str);
                    }
                }
                j++;
                str = strtok(NULL, " ");
            }
            k++;
        }
        fclose(mol_file);
    }
}

void get_differences(struct Job job, struct CharPair mol_nums[job.molecule_types_number],
                     int mol_lengths[job.molecule_types_number],
                     struct Atom differences[job.molecule_types_number][job.max_mol_num][job.max_mol_length],
                     struct Atom molecules[job.molecule_types_number][job.max_mol_length])
{
    float mean[job.molecule_types_number][3];
    int i;
    for(i = 0; i < job.molecule_types_number; i++)
    {
        if (mol_lengths[i] > 1)
        {
            int j;
            for(j = 0; j < 3; j++)
            {
                mean[i][j] = 0.;
            }
            for(j = 0; j < mol_lengths[i]; j++)
            {
                mean[i][0] += molecules[i][j].x;
                mean[i][1] += molecules[i][j].y;
                mean[i][2] += molecules[i][j].z;
            }
            for(j = 0; j < 3; j++)
            {
                mean[i][j] /= mol_lengths[i];
            }
        }
    }
    for(i = 0; i < job.molecule_types_number; i++)
    {
        if (mol_lengths[i] > 1)
        {
            int j;
            for(j = 0; j < atoi(mol_nums[i].keyword); j++)
            {
                int k;
                for(k = 0; k < mol_lengths[i]; k++)
                {
                    strncpy(differences[i][j][k].label, molecules[i][k].label, 3);
                    differences[i][j][k].x = molecules[i][k].x - mean[i][0];
                    differences[i][j][k].y = molecules[i][k].y - mean[i][1];
                    differences[i][j][k].z = molecules[i][k].z - mean[i][2];
                    int r = 0;
                    for(r = 0; r < job.scattering_data_number; r++)
                    {
                        differences[i][j][k].sl[r] = molecules[i][k].sl[r];
                    }
                }
            }
        }
        else
        {
            int j;
            for(j = 0; j < atoi(mol_nums[i].keyword); j++)
            {
                int k;
                for(k = 0; k < mol_lengths[i]; k++)
                {
                    strncpy(differences[i][j][k].label, molecules[i][k].label, 3);
                    differences[i][j][k].x = 0.;
                    differences[i][j][k].y = 0.;
                    differences[i][j][k].z = 0.;
                    int r = 0;
                    for(r = 0; r < job.scattering_data_number; r++)
                    {
                        differences[i][j][k].sl[r] = molecules[i][k].sl[r];
                    }
                }
            }
        }
    }
}

int get_pop_per_core(int population, int rank)
{
    return population / rank;
}

float rand_float(float max, float min)
{
    float x = (float)rand() / ((float)RAND_MAX / max) + min;
    return x;
}


void build_population(struct Job job,
                      struct PosAng population[job.population_per_core][job.molecule_types_number][job.max_mol_num],
                      struct CharPair mol_nums[job.molecule_types_number])
{
    int i;
    for (i = 0; i < job.population_per_core; i++)
    {
        int j;
        for (j = 0; j < job.molecule_types_number; j++)
        {
            int k;
            for (k = 0; k < atoi(mol_nums[j].keyword); k++)
            {
                int l;
                for (l = 0; l < 3; l++)
                {
                    population[i][j][k].position[0] = rand_float(job.cell[0], 0.);
                    population[i][j][k].position[1] = rand_float(job.cell[1], 0.);
                    population[i][j][k].position[2] = rand_float(job.cell[2], 0.);
                    population[i][j][k].angle[0] = rand_float(2 * PI, 0.);
                    population[i][j][k].angle[1] = rand_float(2 * PI, 0.);
                    population[i][j][k].angle[2] = rand_float(2 * PI, 0.);
                }
            }
        }
    }
}


void rotation_matrix(float x, float y, float z, float rm[3][3])
{
    rm[0][0] = cos(x) * cos(y) * cos(z) - sin(x) * sin(z);
    rm[0][1] = -cos(x) * cos(y) * sin(z) - sin(x) * cos(z);
    rm[0][2] = cos(x) * sin(y);
    rm[1][0] = sin(x) * cos(y) * cos(z) + cos(x) * sin(z);
    rm[1][1] = -sin(x) * cos(y) * sin(z) + cos(x) * cos(z);
    rm[1][2] = sin(x) * sin(y);
    rm[2][0] = -sin(y) * cos(z);
    rm[2][1] = sin(y) * sin(z);
    rm[2][2] = cos(y);
}

void convert_to_atomic(struct Job job, struct CharPair mol_nums[job.molecule_types_number],
                       struct PosAng population[job.population_per_core][job.molecule_types_number][job.max_mol_num],
                       int mol_lengths[job.molecule_types_number],
                       struct Atom atomic[job.population_per_core][job.molecule_types_number][job.max_mol_num][job.max_mol_length],
                       struct Atom differences[job.molecule_types_number][job.max_mol_num][job.max_mol_length])
{
    int i;
    for (i = 0; i < job.population_per_core; i++)
    {
        int j;
        for (j = 0; j < job.molecule_types_number; j++)
        {
            int k;
            for (k = 0; k < atoi(mol_nums[j].keyword); k++)
            {
                float rm[3][3];
                rotation_matrix(population[i][j][k].angle[0], population[i][j][k].angle[1],
                                population[i][j][k].angle[2], rm);
                int l;
                for (l = 0; l < mol_lengths[j]; l++)
                {
                    strncpy(atomic[i][j][k][l].label, differences[j][k][l].label, 3);
                    atomic[i][j][k][l].x = population[i][j][k].position[0] - differences[j][k][l].x;
                    atomic[i][j][k][l].y = population[i][j][k].position[1] - differences[j][k][l].y;
                    atomic[i][j][k][l].z = population[i][j][k].position[2] - differences[j][k][l].z;
                    int r = 0;
                    for (r = 0; r < job.scattering_data_number; r++) {
                        atomic[i][j][k][l].sl[r] = differences[j][k][l].sl[r];
                    }
                }
                struct Atom other[mol_lengths[j]];
                for (l = 0; l < mol_lengths[j]; l++)
                {
                    other[l].x = atomic[i][j][k][l].x - population[i][j][k].position[0];
                    other[l].y = atomic[i][j][k][l].y - population[i][j][k].position[1];
                    other[l].z = atomic[i][j][k][l].z - population[i][j][k].position[2];
                }
                struct Atom another[mol_lengths[j]];
                for (l = 0; l < mol_lengths[j]; l++)
                {
                    another[l].x = rm[0][0] * other[l].x + rm[0][1] * other[l].y + rm[0][2] * other[l].z;
                    another[l].y = rm[1][0] * other[l].x + rm[1][1] * other[l].y + rm[1][2] * other[l].z;
                    another[l].z = rm[2][0] * other[l].x + rm[2][1] * other[l].y + rm[2][2] * other[l].z;
                }
                for (l = 0; l < mol_lengths[j]; l++)
                {
                    atomic[i][j][k][l].x = another[l].x + population[i][j][k].position[0];
                    atomic[i][j][k][l].y = another[l].y + population[i][j][k].position[1];
                    atomic[i][j][k][l].z = another[l].z + population[i][j][k].position[2];
                }
            }
        }
    }
}


void write_to_xyz(struct Job job,
                  struct PosAng population[job.population_per_core][job.molecule_types_number][job.max_mol_num],
                  struct Atom differences[job.molecule_types_number][job.max_mol_num][job.max_mol_length],
                  struct CharPair mol_nums[job.molecule_types_number], int mol_lengths[job.molecule_types_number],
                  int to_print, int ranking, int iter)
{
    struct Atom atomic[job.population_per_core][job.molecule_types_number][job.max_mol_num][job.max_mol_length];
    convert_to_atomic(job, mol_nums, population, mol_lengths, atomic, differences);
    char str[50];
    sprintf(str, "fitoogbest_%d_%d.xyz", iter, ranking);
    FILE *f = fopen(str, "w");
    if (f == NULL)
    {
        printf("Error opening output xyz file!\n");
    }
    else
    {
        int num_atoms = 0;
        int i;
        for (i = 0; i < job.molecule_types_number; i++)
        {
            int j;
            for (j = 0; j < atoi(mol_nums[i].keyword); j++)
            {
                int k;
                for (k = 0; k < mol_lengths[i]; k++)
                {
                    num_atoms += 1;
                }
            }
        }
        fprintf(f, "%d\nThis xyz file was produced by fitoog.\n", num_atoms);
        for (i = 0; i < job.molecule_types_number; i++)
        {
            int j;
            for (j = 0; j < atoi(mol_nums[i].keyword); j++)
            {
                int k;
                for (k = 0; k < mol_lengths[i]; k++)
                {
                    fprintf(f, "%s %f %f %f \n", atomic[to_print][i][j][k].label, atomic[to_print][i][j][k].x,
                            atomic[to_print][i][j][k].y, atomic[to_print][i][j][k].z);
                }
            }
        }
        fclose(f);
    }
}

void normalise(struct Job job, struct Data exp_data[job.scattering_data_number][job.max_data_length],
               struct Data sim_data[job.scattering_data_number][job.max_data_length],
               int data_lengths[job.scattering_data_number], int i)
{
    int random = rand() % data_lengths[i];
    float norm = exp_data[i][random].i / sim_data[i][random].i;
    int j;
    for (j = 0; j < data_lengths[i]; j++)
    {
        sim_data[i][j].i = sim_data[i][j].i * norm;
    }
}

void diffraction_calculator(struct Job job,
                            struct Atom atomic[job.molecule_types_number][job.max_mol_num][job.max_mol_length],
                            struct Data exp_data[job.scattering_data_number][job.max_data_length],
                            struct Data sim_data[job.scattering_data_number][job.max_data_length],
                            int data_lengths[job.scattering_data_number],
                            struct CharPair mol_nums[job.molecule_types_number],
                            int mol_lengths[job.molecule_types_number])
{
    int i;
    for (i = 0; i < job.scattering_data_number; i++)
    {
        int q;
        for (q = 0; q < data_lengths[i]; q++)
        {
            float sum_int = 0.;
            float k_min = (-1. * ((float) job.golden_vectors - 1.)) / 2.;
            float k_max = ((float) job.golden_vectors - 1.) / 2.;
            float a;
            for (a = k_min; a < (k_max + 1); a += 1)
            {
                float qx = exp_data[i][q].q * cos(asin((2. * a) / (float) job.golden_vectors)) * cos((2. * PI * a) / (PHI));
                float qy = exp_data[i][q].q * cos(asin((2. * a) / (float) job.golden_vectors)) * sin((2. * PI * a) / (PHI));
                float qz = (2. * a * exp_data[i][q].q) / (float) job.golden_vectors;
                float sum_coss = 0.;
                float sum_sinn = 0.;
                int j;
                for (j = 0; j < job.molecule_types_number; j++)
                {
                    int k;
                    for (k = 0; k < atoi(mol_nums[j].keyword); k++)
                    {
                        int l;
                        for (l = 0; l < mol_lengths[j]; l++)
                        {
                            float dot = qx * atomic[j][k][l].x + qy * atomic[j][k][l].y + qz * atomic[j][k][l].z;
                            float coss = atomic[j][k][l].sl[i] * cos(dot);
                            float sinn = atomic[j][k][l].sl[i] * sin(dot);
                            sum_coss += coss;
                            sum_sinn += sinn;
                        }
                        float sq_coss = sum_coss * sum_coss;
                        float sq_sinn = sum_sinn * sum_sinn;
                        sum_int += (sq_coss + sq_sinn);
                    }
                }
            }
            float inten = sum_int / (float) job.golden_vectors;
            sim_data[i][q].i = inten;
            sim_data[i][q].q = exp_data[i][q].q;
        }
        normalise(job, exp_data, sim_data, data_lengths, i);
    }
}

float min(int length, float array[length])
{
    float a = array[0];
    int i;
    for (i = 0; i < length; i++)
    {
        if (array[i] < a)
        {
            a = array[i];
        }
    }
    return a;
}

float get_chisq(struct Job job, struct Data exp_data[job.scattering_data_number][job.max_data_length],
               struct Data sim_data[job.scattering_data_number][job.max_data_length],
               int data_lengths[job.scattering_data_number])
{
    int j;
    float chi_sq = 0.;
    float min_cell = min(3, job.cell);
    for (j = 0; j < job.scattering_data_number; j++)
    {
        int k;
        for (k = 0; k < data_lengths[j]; k++)
        {
            if (exp_data[j][k].q >= (2 * PI) / (min_cell / 2.))
            {
                chi_sq += (exp_data[j][k].i - sim_data[j][k].i) *
                        (exp_data[j][k].i - sim_data[j][k].i) / exp_data[j][k].di;
            }
        }
    }
    return chi_sq;
}


void analyse(struct Job job,
             struct PosAng population[job.population_per_core][job.molecule_types_number][job.max_mol_num],
             struct CharPair mol_nums[job.molecule_types_number], int mol_lengths[job.molecule_types_number],
             struct Atom differences[job.molecule_types_number][job.max_mol_num][job.max_mol_length],
             struct Data exp_data[job.scattering_data_number][job.max_data_length],
             struct Data sim_data[job.scattering_data_number][job.max_data_length],
             int data_lengths[job.scattering_data_number], float chi_sq[job.population_per_core])
{
    struct Atom atomic[job.population_per_core][job.molecule_types_number][job.max_mol_num][job.max_mol_length];
    convert_to_atomic(job, mol_nums, population, mol_lengths, atomic, differences);
    int i;
    for (i = 0; i < job.population_per_core; i++)
    {
        diffraction_calculator(job, atomic[i], exp_data, sim_data, data_lengths, mol_nums, mol_lengths);
        chi_sq[i] = get_chisq(job, exp_data, sim_data, data_lengths);
    }
}


MPI_Comm mpstart(int *n_procs, int *rank)
{
    MPI_Init(NULL, NULL);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, n_procs);
    MPI_Comm_rank(comm, rank);
    return comm;
}


void best_chi_sq_local(struct Job job, float all_chi_sq[job.population_per_core], int best[2], int i)
{
    float a = all_chi_sq[0];
    int besti = i;
    int bestj = 0;
    int j;
    for (j = 0; j < job.population_per_core; j++)
    {
        if (all_chi_sq[j] < a)
        {
            a = all_chi_sq[j];
            bestj = j;
        }
    }
    best[0] = besti;
    best[1] = bestj;
}

void best_chi_sq(struct Job job, int n_procs, float all_chi_sq[n_procs][job.population_per_core], int best[2],
                 double *gbest_chisq)
{
    int besti = 0;
    int bestj = 0;
    int i;
    for (i = 0; i < n_procs; i++)
    {
        int j;
        for (j = 0; j < job.population_per_core; j++)
        {
            if (all_chi_sq[i][j] < *gbest_chisq)
            {
                *gbest_chisq = all_chi_sq[i][j];
                besti = i;
                bestj = j;
            }
        }
    }
    best[0] = besti;
    best[1] = bestj;
}

void initialise_velocities(struct Job job,
                           struct PosAng velocity[job.population_per_core][job.molecule_types_number][job.max_mol_num],
                           int mol_lengths[job.molecule_types_number])
{
    int i;
    for (i = 0; i < job.population_per_core; i++)
    {
        int j;
        for (j = 0; j < job.molecule_types_number; j++)
        {
            int k;
            for (k = 0; k < mol_lengths[j]; k++)
            {
                int l;
                for (l = 0; l < 3; l++) {
                    velocity[i][j][k].position[l] = rand_float(1., 0.);
                    velocity[i][j][k].angle[l] = rand_float(1., 0.);
                }
            }
        }
    }
}

void integrator(struct Job job, int n_procs, int rank,
                   struct PosAng population[job.population_per_core][job.molecule_types_number][job.max_mol_num],
                   struct PosAng pbest[job.population_per_core][job.molecule_types_number][job.max_mol_num],
                   struct PosAng gbest[job.molecule_types_number][job.max_mol_num],
                   int mol_lengths[job.molecule_types_number],
                   struct PosAng velocity[job.population_per_core][job.molecule_types_number][job.max_mol_num]) {
    int i;
    for (i = 0; i < job.population_per_core; i++) {
        int j;
        for (j = 0; j < job.molecule_types_number; j++) {
            int k;
            for (k = 0; k < mol_lengths[j]; k++) {
                int l;
                for (l = 0; l < 3; l++) {
                    velocity[i][j][k].position[l] = velocity[i][j][k].position[l] +
                                                    (rand_float(job.phi[0], 0.0) * (pbest[i][j][k].position[l] -
                                                                                population[i][j][k].position[l])) +
                                                    (rand_float(job.phi[1], 0.0) * (gbest[j][k].position[l] -
                                                                              population[i][j][k].position[l]));
                    velocity[i][j][k].angle[l] = velocity[i][j][k].angle[l] +
                                                 (rand_float(job.phi[0], 0.0) * (pbest[i][j][k].angle[l] -
                                                                             population[i][j][k].angle[l])) +
                                                 (rand_float(job.phi[1], 0.0) * (gbest[j][k].angle[l] -
                                                                           population[i][j][k].angle[l]));
                }
            }
        }
    }
    for (i = 0; i < job.population_per_core; i++) {
        int j;
        for (j = 0; j < job.molecule_types_number; j++) {
            int k;
            for (k = 0; k < mol_lengths[j]; k++) {
                int l;
                for (l = 0; l < 3; l++) {
                    population[i][j][k].position[l] = population[i][j][k].position[l] + velocity[i][j][k].position[l];
                    population[i][j][k].position[l] = fabs(fmod(population[i][j][k].position[l], job.cell[l]));
                    population[i][j][k].angle[l] = population[i][j][k].angle[l] + velocity[i][j][k].angle[l];
                    population[i][j][k].angle[l] = fabs(fmod(population[i][j][k].angle[l], (2. * PI)));
                }
            }
        }
    }
}

void initialise_best(struct Job job,
                     struct PosAng pbest[job.population_per_core][job.molecule_types_number][job.max_mol_num],
                     struct PosAng gbest[job.molecule_types_number][job.max_mol_num],
                     struct PosAng population[job.population_per_core][job.molecule_types_number][job.max_mol_num],
                     int mol_lengths[job.molecule_types_number], double pbest_chisq[job.population_per_core],
                     double *gbest_chisq)
{
    int i;
    for (i = 0; i < job.population_per_core; i++)
    {
        pbest_chisq[i] = 1E20;
        int j;
        for (j = 0; j < job.molecule_types_number; j++)
        {
            int k;
            for (k = 0; k < mol_lengths[j]; k++)
            {
                int l;
                for (l = 0; l < 3; l++) {
                    pbest[i][j][k].position[l] = population[i][j][k].position[l];
                    pbest[i][j][k].angle[l] = population[i][j][k].angle[l];
                }
            }
        }
    }
    *gbest_chisq = 1E20;
    int j;
    for (j = 0; j < job.molecule_types_number; j++)
    {
        int k;
        for (k = 0; k < mol_lengths[j]; k++)
        {
            int l;
            for (l = 0; l < 3; l++) {
                gbest[j][k].position[l] = 0.;
                gbest[j][k].angle[l] = 0.;
            }
        }
    }
}

void update_pbest(struct Job job,
                  struct PosAng pbest[job.population_per_core][job.molecule_types_number][job.max_mol_num],
                  struct PosAng population[job.population_per_core][job.molecule_types_number][job.max_mol_num],
                  int mol_lengths[job.molecule_types_number], double pbest_chisq[job.population_per_core],
                  float chi_sq[job.population_per_core])
{
    int i;
    for (i = 0; i < job.population_per_core; i++)
    {
        if (chi_sq[i] < pbest_chisq[i]) {
            pbest_chisq[i] = chi_sq[i];
            int j;
            for (j = 0; j < job.molecule_types_number; j++) {
                int k;
                for (k = 0; k < mol_lengths[j]; k++) {
                    int l;
                    for (l = 0; l < 3; l++) {
                        pbest[i][j][k].position[l] = population[i][j][k].position[l];
                        pbest[i][j][k].angle[l] = population[i][j][k].angle[l];
                    }
                }
            }
        }
    }
}

void update_gbest(struct Job job, struct PosAng gbest[job.molecule_types_number][job.max_mol_num],
                  struct PosAng population[job.population_per_core][job.molecule_types_number][job.max_mol_num],
                  int mol_lengths[job.molecule_types_number], double *gbest_chisq, int n_procs,
                  float all_chi_sq[n_procs][job.population_per_core], int rank, MPI_Comm comm,
                  struct Atom differences[job.molecule_types_number][job.max_mol_num][job.max_mol_length],
                  struct CharPair mol_nums[job.molecule_types_number], int iter)
{
    int best[2];
    float storage;
    float storage2;
    if (rank == 0){
        best_chi_sq(job, n_procs, all_chi_sq, best, gbest_chisq);
        printf("%d %d %f %f\n", best[0], best[1], all_chi_sq[best[0]][best[1]], *gbest_chisq);
        storage = *gbest_chisq;
        storage2 = all_chi_sq[best[0]][best[1]];
    }
    MPI_Bcast(&best, 2, MPI_INT, 0, comm);
    MPI_Bcast(&storage, 1, MPI_FLOAT, 0, comm);
    MPI_Bcast(&storage2, 1, MPI_FLOAT, 0, comm);
    float best_of_pop[job.molecule_types_number][job.max_mol_num][6];
    if (rank == best[0])
    {
        int i;
        for (i = 0; i < job.molecule_types_number; i++)
        {
            int j;
            for (j = 0; j < mol_lengths[i]; j++)
            {
                int k;
                for (k = 0; k < 3; k++)
                {
                    best_of_pop[i][j][k] = population[best[1]][i][j].position[k];
                    best_of_pop[i][j][k + 3] = population[best[1]][i][j].angle[k + 3];
                }
            }
        }
    }
    if (storage2 <= storage)
    {
        MPI_Bcast(&best_of_pop, job.molecule_types_number * job.max_mol_num * 6, MPI_FLOAT, best[0], comm);
        int i;
        for (i = 0; i < job.molecule_types_number; i++)
        {
            int j;
            for (j = 0; j < mol_lengths[i]; j++)
            {
                int k;
                for (k = 0; k < 3; k++)
                {
                    gbest[i][j].position[k] = best_of_pop[i][j][k];
                    gbest[i][j].angle[k+3] = best_of_pop[i][j][k + 3];
                }
            }
        }
        if (rank == 0)
        {
            write_to_xyz(job, population, differences, mol_nums, mol_lengths, best[1], 0, iter);
        }
    }
    else
    {
        MPI_Bcast(&best_of_pop, job.molecule_types_number * job.max_mol_num * 6, MPI_FLOAT, best[0], comm);
        int i;
        for (i = 0; i < job.molecule_types_number; i++)
        {
            int j;
            for (j = 0; j < mol_lengths[i]; j++)
            {
                int k;
                for (k = 0; k < 3; k++)
                {
                    gbest[i][j].position[k] = gbest[i][j].position[k];
                    gbest[i][j].angle[k+3] = gbest[i][j].angle[k+3];
                }
            }
        }
    }
}

int main(int argc, char *argv[])
{
    double t1, t2;
    t1 = MPI_Wtime();
    int rank, n_procs;
    MPI_Comm comm = mpstart(&n_procs, &rank);
    srand(rank + 1 * time(NULL));

    struct Job job;
    read_input_file(&job);
    job.population_per_core = get_pop_per_core(job.population_size, n_procs);

    struct CharPair mol_types[job.molecule_types_number];
    get_iterator_labels(job.molecule_types_number, mol_types, "molecule");

    struct CharPair data_types[job.scattering_data_number];
    get_iterator_labels(job.scattering_data_number, data_types, "data");

    struct CharPair mol_nums[job.molecule_types_number];
    get_iterator_labels(job.molecule_types_number, mol_nums, "num_mole");
    job.max_mol_num = find_max_charpair(job.molecule_types_number, mol_nums);

    int mol_lengths[job.molecule_types_number];
    get_files_lengths(job.molecule_types_number, mol_lengths, mol_types);
    job.max_mol_length = find_max_int(job.molecule_types_number, mol_lengths);

    int data_lengths[job.scattering_data_number];
    get_files_lengths(job.scattering_data_number, data_lengths, data_types);
    job.max_data_length = find_max_int(job.scattering_data_number, data_lengths);

    struct Data exp_data[job.scattering_data_number][job.max_data_length];
    struct Data sim_data[job.scattering_data_number][job.max_data_length];
    get_data_points(job, exp_data, data_types);

    struct Atom molecules[job.molecule_types_number][job.max_mol_length];
    struct Atom differences[job.molecule_types_number][job.max_mol_num][job.max_mol_length];
    struct PosAng population[job.population_per_core][job.molecule_types_number][job.max_mol_num];

    if (job.restart == 0)
    {
        get_atom_positions(job, molecules, mol_types);

        get_differences(job, mol_nums, mol_lengths, differences, molecules);

        build_population(job, population, mol_nums);
    }
    else
    {
        exit(1);
    }

    struct PosAng velocity[job.population_per_core][job.molecule_types_number][job.max_mol_num];
    struct PosAng pbest[job.population_per_core][job.molecule_types_number][job.max_mol_num];
    struct PosAng gbest[job.molecule_types_number][job.max_mol_num];
    double pbest_chisq[job.population_per_core];
    double gbest_chisq;

    initialise_best(job, pbest, gbest, population, mol_lengths, pbest_chisq, &gbest_chisq);
    initialise_velocities(job, velocity, mol_lengths);

    int i;
    for (i = 0; i < job.steps_number; i++)
    {
        float chi_sq[job.population_per_core];
        analyse(job, population, mol_nums, mol_lengths, differences, exp_data, sim_data, data_lengths, chi_sq);

        float all_chi_sq[n_procs][job.population_per_core];
        MPI_Gather(chi_sq, job.population_per_core, MPI_FLOAT, all_chi_sq, job.population_per_core, MPI_FLOAT, 0, comm);

        update_pbest(job, pbest, population, mol_lengths, pbest_chisq, chi_sq);
        update_gbest(job, gbest, population, mol_lengths, &gbest_chisq, n_procs, all_chi_sq, rank, comm, differences,
                     mol_nums, i);

        integrator(job, n_procs, rank, population, pbest, gbest, mol_lengths, velocity);
    }

    t2 = MPI_Wtime();
    if (rank == 0)
    {
        printf( "Elapsed time is %f\n", t2 - t1 );
    }
    MPI_Finalize();
}
