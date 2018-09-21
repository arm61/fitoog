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
    int index;
    char label[3];
    float x;
    float y;
    float z;
    float sl[10];
};

struct Bond
{
    int molecule_index;
    int a;
    int b;
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
    FILE *restart_file;
    int em_freq;
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

int read_int(char line[100], const char wipe[50])
{
    char temp[50];
    strncpy(temp, wipe, 50);
    strncpy(temp, first_char_after_space(line), 50);
    strcat(temp, "\0");
    return atoi(temp);
}

float read_float(char line[100], const char wipe[50])
{
    char temp[50];
    strncpy(temp, wipe, 50);
    strncpy(temp, first_char_after_space(line), 50);
    strcat(temp, "\0");
    return atof(temp);
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
            job->molecule_types_number = read_int(line, wipe);
        }
        found = find_word_in_line(line, "phi0");
        if(found != NULL)
        {
            job->phi[0] = read_float(line, wipe);
        }
        found = find_word_in_line(line, "phi1");
        if(found != NULL)
        {
            job->phi[1] = read_float(line, wipe);
        }
        found = find_word_in_line(line, "cell0");
        if(found != NULL)
        {
            job->cell[0] = read_float(line, wipe);
        }
        found = find_word_in_line(line, "cell1");
        if(found != NULL)
        {
            job->cell[1] = read_float(line, wipe);
        }
        found = find_word_in_line(line, "cell2");
        if(found != NULL)
        {
            job->cell[2] = read_float(line, wipe);
        }
        found = find_word_in_line(line, "num_data");
        if(found != NULL)
        {
            job->scattering_data_number = read_int(line, wipe);
        }
        found = find_word_in_line(line, "pop_size");
        if(found != NULL)
        {
            job->population_size = read_int(line, wipe);
        }
        found = find_word_in_line(line, "num_steps");
        if(found != NULL)
        {
            job->steps_number = read_int(line, wipe);
        }
        found = find_word_in_line(line, "restart");
        if(found != NULL)
        {
            if(first_char_after_space(line)[0]=='y')
            {
                job->restart = 1;
            }
            else
            {
                job->restart = 0;
            }
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
                    molecules[i][k].index = atoi(str);
                }
                if (j == 1)
                {
                    strncpy(molecules[i][k].label, str, 3);
                }
                if (j == 2)
                {
                    molecules[i][k].x = atof(str);
                }
                if (j == 3)
                {
                    molecules[i][k].y = atof(str);
                }
                if (j == 4)
                {
                    molecules[i][k].z = atof(str);
                }
                int r;
                for (r = 0; r < job.scattering_data_number; r++)
                {
                    if (j == r + 5)
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
                     struct Atom differences[job.population_per_core][job.molecule_types_number][job.max_mol_num][job.max_mol_length],
                     struct Atom molecules[job.molecule_types_number][job.max_mol_length])
{
    float mean[job.molecule_types_number][3];
    int p;
    for (p = 0; p < job.population_per_core; p++)
    {
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
                        differences[p][i][j][k].index = molecules[i][k].index;
                        strncpy(differences[p][i][j][k].label, molecules[i][k].label, 3);
                        differences[p][i][j][k].x = molecules[i][k].x - mean[i][0];
                        differences[p][i][j][k].y = molecules[i][k].y - mean[i][1];
                        differences[p][i][j][k].z = molecules[i][k].z - mean[i][2];
                        int r = 0;
                        for(r = 0; r < job.scattering_data_number; r++)
                        {
                            differences[p][i][j][k].sl[r] = molecules[i][k].sl[r];
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
                        differences[p][i][j][k].index = molecules[i][k].index;
                        strncpy(differences[p][i][j][k].label, molecules[i][k].label, 3);
                        differences[p][i][j][k].x = 0.;
                        differences[p][i][j][k].y = 0.;
                        differences[p][i][j][k].z = 0.;
                        int r = 0;
                        for(r = 0; r < job.scattering_data_number; r++)
                        {
                            differences[p][i][j][k].sl[r] = molecules[i][k].sl[r];
                        }
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
                       struct Atom differences[job.population_per_core][job.molecule_types_number][job.max_mol_num][job.max_mol_length])
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
                    strncpy(atomic[i][j][k][l].label, differences[i][j][k][l].label, 3);
                    atomic[i][j][k][l].x = population[i][j][k].position[0] - differences[i][j][k][l].x;
                    atomic[i][j][k][l].y = population[i][j][k].position[1] - differences[i][j][k][l].y;
                    atomic[i][j][k][l].z = population[i][j][k].position[2] - differences[i][j][k][l].z;
                    int r = 0;
                    for (r = 0; r < job.scattering_data_number; r++) {
                        atomic[i][j][k][l].sl[r] = differences[i][j][k][l].sl[r];
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

void write_to_fit(struct Job job,
                  struct PosAng population[job.population_per_core][job.molecule_types_number][job.max_mol_num],
                  struct Atom differences[job.population_per_core][job.molecule_types_number][job.max_mol_num][job.max_mol_length],
                  struct CharPair mol_nums[job.molecule_types_number], int mol_lengths[job.molecule_types_number],
                  int to_print, int ranking, int iter)
{
    struct Atom atomic[job.population_per_core][job.molecule_types_number][job.max_mol_num][job.max_mol_length];
    convert_to_atomic(job, mol_nums, population, mol_lengths, atomic, differences);
    char str[50];
    sprintf(str, "fitoogbest_%d_%d.fit", iter, ranking);
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
        int count = 0;
        for (i = 0; i < job.molecule_types_number; i++)
        {
            int j;
            for (j = 0; j < atoi(mol_nums[i].keyword); j++)
            {
                int k;
                for (k = 0; k < mol_lengths[i]; k++)
                {
                    fprintf(f, "%d %d %s %f %f %f ", count, k, atomic[to_print][i][j][k].label, atomic[to_print][i][j][k].x,
                            atomic[to_print][i][j][k].y, atomic[to_print][i][j][k].z);
                    int r;
                    for (r = 0; r < job.scattering_data_number; r++)
                    {
                        fprintf(f, "%f ", atomic[to_print][i][j][k].sl[r]);
                    }
                    fprintf(f, "\n");
                }
                count += 1;
            }
        }
        fclose(f);
    }
}

void write_to_xyz(struct Job job,
                  struct PosAng population[job.population_per_core][job.molecule_types_number][job.max_mol_num],
                  struct Atom differences[job.population_per_core][job.molecule_types_number][job.max_mol_num][job.max_mol_length],
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
             struct Atom differences[job.population_per_core][job.molecule_types_number][job.max_mol_num][job.max_mol_length],
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

void get_energy(struct Job job, struct Atom atomic[job.population_per_core][job.molecule_types_number][job.max_mol_num][job.max_mol_length], int to_minim, struct CharPair mol_nums[job.molecule_types_number], int mol_lengths[job.molecule_types_number], struct Bond bonds[job.max_mol_length], int num_bonds, int num_atoms, float force_x[num_atoms], float force_y[num_atoms], float force_z[num_atoms], float *total_energy, float *total_force)
{
    int i;
    int pairs = 0;

    float sigma6 = pow(6.2, 6);
    float sigma12 = sigma6 * sigma6;

    // intramolecular interactions
    for (i = 0; i < job.molecule_types_number; i++)
    {
        int j;
        for (j = 0; j < atoi(mol_nums[i].keyword); j++)
        {
            int k;
            for (k = 0; k < mol_lengths[i] - 1; k++)
            {
                int k2;
                for (k2 = k + 1; k2 < mol_lengths[i]; k2++)
                {
                    pairs += 1;
                    int kk = (i*atoi(mol_nums[i].keyword)*mol_lengths[i-1])+(j*mol_lengths[i])+k;
                    int kk2 = (i*atoi(mol_nums[i].keyword)*mol_lengths[i-1])+(j*mol_lengths[i])+k2;
                    float dx = atomic[to_minim][i][j][k2].x - atomic[to_minim][i][j][k].x;
                    float dy = atomic[to_minim][i][j][k2].y - atomic[to_minim][i][j][k].y;
                    float dz = atomic[to_minim][i][j][k2].z - atomic[to_minim][i][j][k].z;
                    if (fabs(dx) > 0.5 * job.cell[0])
                    {
                        dx *= 1 - job.cell[0] / fabs(dx);
                    }
                    if (fabs(dy) > 0.5 * job.cell[1])
                    {
                        dy *= 1 - job.cell[1] / fabs(dy);
                    }
                    if (fabs(dz) > 0.5 * job.cell[2])
                    {
                        dz *= 1 - job.cell[2] / fabs(dz);
                    }
                    float d = sqrt(dx * dx + dy * dy + dz * dz);
                    if (d <= 15.)
                    {
                        float d2 = d * d;
                        float d6 = d2 * d2 * d2;
                        float d12 = d6 * d6;
                        float energy = 4 * 0.000332 * (sigma12 / d12 - sigma6 / d6);
                        float force = 4 * 0.000332 * ((6 * sigma6) / (d6 * d) - (12 * sigma12) / (d12 * d));
                        *total_energy += energy;
                        *total_force += force;
                        force_x[kk] += force * dx / d;
                        force_x[kk2] -= force * dx / d;
                        force_y[kk] += force * dy / d;
                        force_y[kk2] -= force * dy / d;
                        force_z[kk] += force * dz / d;
                        force_z[kk2] -= force * dz / d;
                    }
                    int l;
                    for (l = 0; l < num_bonds; l++)
                    {
                        if (bonds[l].molecule_index == i && bonds[l].a == k && bonds[l].b == k2){
                            float energy = 0.0020757 * 0.5 * pow((d - 4.7), 2);
                            float force = 0.0020757  * (d - 4.7);
                            *total_energy += energy;
                            *total_force += force;
                            force_x[kk] += force * dx / d;
                            force_x[kk2] -= force * dx / d;
                            force_y[kk] += force * dy / d;
                            force_y[kk2] -= force * dy / d;
                            force_z[kk] += force * dz / d;
                            force_z[kk2] -= force * dz / d;
                        }
                    }
                }
            }
        }
    }
    // intermolecular interactions for same molecule types
    for (i = 0; i < job.molecule_types_number; i++)
    {
        int j;
        for (j = 0; j < atoi(mol_nums[i].keyword) - 1; j++)
        {
            int j2;
            for (j2 = j + 1; j2 <  atoi(mol_nums[i].keyword); j2++)
            {
                int k;
                for (k = 0; k < mol_lengths[i]; k++)
                {
                    int k2;
                    for (k2 = 0; k2 < mol_lengths[i]; k2++)
                    {
                        pairs += 1;
                        int kk = (i*atoi(mol_nums[i].keyword)*mol_lengths[i-1])+(j*mol_lengths[i])+k;
                        int kk2 = (i*atoi(mol_nums[i].keyword)*mol_lengths[i-1])+(j2*mol_lengths[i])+k2;
                        float dx = atomic[to_minim][i][j2][k2].x - atomic[to_minim][i][j][k].x;
                        float dy = atomic[to_minim][i][j2][k2].y - atomic[to_minim][i][j][k].y;
                        float dz = atomic[to_minim][i][j2][k2].z - atomic[to_minim][i][j][k].z;
                        if (fabs(dx) > 0.5 * job.cell[0])
                        {
                            dx *= 1 - job.cell[0] / fabs(dx);
                        }
                        if (fabs(dy) > 0.5 * job.cell[1])
                        {
                            dy *= 1 - job.cell[1] / fabs(dy);
                        }
                        if (fabs(dz) > 0.5 * job.cell[2])
                        {
                            dz *= 1 - job.cell[2] / fabs(dz);
                        }
                        float d = sqrt(dx * dx + dy * dy + dz * dz);
                        if (d <= 15.)
                        {
                            float d2 = d * d;
                            float d6 = d2 * d2 * d2;
                            float d12 = d6 * d6;
                            float energy = 4 * 0.000332 * (sigma12 / d12 - sigma6 / d6);
                            float force = 4 * 0.000332 * ((6 * sigma6) / (d6 * d) - (12 * sigma12) / (d12 * d));
                            *total_energy += energy;
                            *total_force += force;
                            force_x[kk] += force * dx / d;
                            force_x[kk2] -= force * dx / d;
                            force_y[kk] += force * dy / d;
                            force_y[kk2] -= force * dy / d;
                            force_z[kk] += force * dz / d;
                            force_z[kk2] -= force * dz / d;
                        }
                    }
                }
            }
        }
    }
    // intermolecular interactions for different molecule types
    for (i = 0; i < job.molecule_types_number - 1; i++)
    {
        int i2;
        for (i2 = i + 1; i2 < job.molecule_types_number; i2++)
        {
            int j;
            for (j = 0; j <  atoi(mol_nums[i].keyword); j++)
            {
                int j2;
                for (j2 = 0; j2 <  atoi(mol_nums[i2].keyword); j2++)
                {
                    int k;
                    for (k = 0; k < mol_lengths[i]; k++)
                    {
                        int k2;
                        for (k2 = 0; k2 < mol_lengths[i2]; k2++)
                        {
                            pairs += 1;
                            int kk = (i*atoi(mol_nums[i].keyword)*mol_lengths[i-1])+(j*mol_lengths[i])+k;
                            int kk2 = (i2*atoi(mol_nums[i2].keyword)*mol_lengths[i2-1])+(j2*mol_lengths[i])+k2;
                            float dx = atomic[to_minim][i2][j2][k2].x - atomic[to_minim][i][j][k].x;
                            float dy = atomic[to_minim][i2][j2][k2].y - atomic[to_minim][i][j][k].y;
                            float dz = atomic[to_minim][i2][j2][k2].z - atomic[to_minim][i][j][k].z;
                            if (fabs(dx) > 0.5 * job.cell[0])
                            {
                                dx *= 1 - job.cell[0] / fabs(dx);
                            }
                            if (fabs(dy) > 0.5 * job.cell[1])
                            {
                                dy *= 1 - job.cell[1] / fabs(dy);
                            }
                            if (fabs(dz) > 0.5 * job.cell[2])
                            {
                                dz *= 1 - job.cell[2] / fabs(dz);
                            }
                            float d = sqrt(dx * dx + dy * dy + dz * dz);
                            if (d <= 15.)
                            {
                                float d2 = d * d;
                                float d6 = d2 * d2 * d2;
                                float d12 = d6 * d6;
                                float energy = 4 * 0.000332 * (sigma12 / d12 - sigma6 / d6);
                                float force = 4 * 0.000332 * ((6 * sigma6) / (d6 * d) - (12 * sigma12) / (d12 * d));
                                *total_energy += energy;
                                *total_force += force;
                                force_x[kk] += force * dx / d;
                                force_x[kk2] -= force * dx / d;
                                force_y[kk] += force * dy / d;
                                force_y[kk2] -= force * dy / d;
                                force_z[kk] += force * dz / d;
                                force_z[kk2] -= force * dz / d;
                            }
                        }
                    }
                }
            }
        }
    }
}

void move(struct Job job, struct Atom atomic[job.population_per_core][job.molecule_types_number][job.max_mol_num][job.max_mol_length], int to_minim, struct CharPair mol_nums[job.molecule_types_number], int mol_lengths[job.molecule_types_number], int num_atoms, float force_x[num_atoms], float force_y[num_atoms], float force_z[num_atoms], float hn)
{
    int i;
    for (i = 0; i < job.molecule_types_number; i++)
    {
        int j;
        for (j = 0; j < atoi(mol_nums[i].keyword); j++)
        {
            int k;
            for (k = 0; k < mol_lengths[i] - 1; k++)
            {
                int kk = (i*atoi(mol_nums[i].keyword)*mol_lengths[i-1])+(j*mol_lengths[i])+k;
                float largest_f = sqrt(force_x[kk] * force_x[kk]);
                if (largest_f < sqrt(force_y[kk] * force_y[kk]))
                {
                    largest_f = sqrt(force_y[kk] * force_y[kk]);
                }
                if (largest_f < sqrt(force_z[kk] * force_z[kk]))
                {
                    largest_f = sqrt(force_z[kk] * force_z[kk]);
                }
                atomic[to_minim][i][j][k].x += force_x[kk] / largest_f * hn;
                atomic[to_minim][i][j][k].y += force_y[kk] / largest_f * hn;
                atomic[to_minim][i][j][k].z += force_z[kk] / largest_f * hn;
            }
        }
    }
}

void get_differences_restart(struct Job job, struct CharPair mol_nums[job.molecule_types_number],
                             int mol_lengths[job.molecule_types_number],
                             struct Atom differences[job.population_per_core][job.molecule_types_number][job.max_mol_num][job.max_mol_length], int to_minim,
                             struct Atom atomic[job.population_per_core][job.molecule_types_number][job.max_mol_num][job.max_mol_length])
{
    float mean[job.molecule_types_number][job.max_mol_num][3];
    int i;
    for(i = 0; i < job.molecule_types_number; i++)
    {
        int j;
        for (j = 0; j < atoi(mol_nums[i].keyword); j++)
        {
            if (mol_lengths[i] > 1)
            {
                int k;
                for(k = 0; k < 3; k++)
                {
                    mean[i][j][k] = 0.;
                }
                for(k = 0; k < mol_lengths[i]; k++)
                {
                    mean[i][j][0] += atomic[to_minim][i][j][k].x;
                    mean[i][j][1] += atomic[to_minim][i][j][k].y;
                    mean[i][j][2] += atomic[to_minim][i][j][k].z;
                }
                for(k = 0; k < 3; k++)
                {
                    mean[i][j][k] /= mol_lengths[i];
                }
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
                    strncpy(differences[to_minim][i][j][k].label, atomic[to_minim][i][j][k].label, 3);
                    differences[to_minim][i][j][k].x = atomic[to_minim][i][j][k].x - mean[i][j][0];
                    differences[to_minim][i][j][k].y = atomic[to_minim][i][j][k].y - mean[i][j][1];
                    differences[to_minim][i][j][k].z = atomic[to_minim][i][j][k].z - mean[i][j][2];
                    int r = 0;
                    for(r = 0; r < job.scattering_data_number; r++)
                    {
                        differences[to_minim][i][j][k].sl[r] = atomic[to_minim][i][j][k].sl[r];
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
                    strncpy(differences[to_minim][i][j][k].label, atomic[to_minim][i][j][k].label, 3);
                    differences[to_minim][i][j][k].x = 0.;
                    differences[to_minim][i][j][k].y = 0.;
                    differences[to_minim][i][j][k].z = 0.;
                    int r = 0;
                    for(r = 0; r < job.scattering_data_number; r++)
                    {
                        differences[to_minim][i][j][k].sl[r] = atomic[to_minim][i][j][k].sl[r];
                    }
                }
            }
        }
    }
}

void update_pop(struct Job job, struct PosAng population[job.population_per_core][job.molecule_types_number][job.max_mol_num],
                struct CharPair mol_nums[job.molecule_types_number], int mol_lengths[job.molecule_types_number], struct Atom atomic[job.population_per_core][job.molecule_types_number][job.max_mol_num][job.max_mol_length])
{
    int i;
    for (i = 0; i < job.population_per_core; i++)
    {
      int j;
      for(j = 0; j < job.molecule_types_number; j++)
      {
          int k;
          for (k = 0; k < atoi(mol_nums[j].keyword); k++)
          {
              if (mol_lengths[j] > 1)
              {
                  int l;
                  for(l = 0; l < 3; l++)
                  {
                      population[i][j][k].position[l] = 0.;
                  }
                  for(l = 0; l < mol_lengths[j]; l++)
                  {
                      population[i][j][k].position[0] += atomic[i][j][k][l].x;
                      population[i][j][k].position[1] += atomic[i][j][k][l].y;
                      population[i][j][k].position[2] += atomic[i][j][k][l].z;
                  }
                  for(l = 0; l < 3; l++)
                  {
                      population[i][j][k].position[l] /= mol_lengths[j];
                  }
              }
              else
              {
                population[i][j][k].position[0] = atomic[i][j][k][0].x;
                population[i][j][k].position[1] = atomic[i][j][k][0].y;
                population[i][j][k].position[2] = atomic[i][j][k][0].z;
              }
          }
      }
    }
    for (i = 0; i < job.population_per_core; i++)
    {
      int j;
      for(j = 0; j < job.molecule_types_number; j++)
      {
          int k;
          for (k = 0; k < atoi(mol_nums[j].keyword); k++)
          {
            if (mol_lengths[j] > 1)
            {
              int l;
                for(l = 0; l < mol_lengths[j]; l++)
                {
                  atomic[i][j][k][l].x -= population[i][j][k].position[0];
                  atomic[i][j][k][l].y -= population[i][j][k].position[1];
                  atomic[i][j][k][l].z -= population[i][j][k].position[2];
                }
            }
            else
            {
              population[i][j][k].angle[0] = 0.;
              population[i][j][k].angle[1] = 0.;
              population[i][j][k].angle[2] = 0.;
            }
          }
        }
      }
}

void energy_minimisation(struct Job job,
                         struct PosAng population[job.population_per_core][job.molecule_types_number][job.max_mol_num],
                         struct Atom differences[job.population_per_core][job.molecule_types_number][job.max_mol_num][job.max_mol_length],
                         struct CharPair mol_nums[job.molecule_types_number], int mol_lengths[job.molecule_types_number], struct Bond bonds[job.max_mol_length], int num_bonds, struct Data exp_data[job.scattering_data_number][job.max_data_length], struct Data sim_data[job.scattering_data_number][job.max_data_length], int data_lengths[job.scattering_data_number])
{
    struct Atom atomic[job.population_per_core][job.molecule_types_number][job.max_mol_num][job.max_mol_length];
    struct Atom old_atomic[job.population_per_core][job.molecule_types_number][job.max_mol_num][job.max_mol_length];

    convert_to_atomic(job, mol_nums, population, mol_lengths, atomic, differences);
    int p;
    for (p = 0; p < job.population_per_core; p++)
    {
        int num_atoms = 0;
        int i;
        for (i = 0; i < job.molecule_types_number; i++)
        {
            num_atoms += atoi(mol_nums[i].keyword) * mol_lengths[i];
        }
        int count;
        float old_force = 0;
        float old_energy = 100000;
        float total_energy = 0;

        double hn = 0.001;
        float total_force = 1000000;
        int count_bad = 0;
        while (abs(total_force) > 1000 && count_bad < 100)
        {
            old_force = total_force;
            old_energy = total_energy;
            total_force = 0;
            total_energy = 0;
            float force_x[num_atoms], force_y[num_atoms], force_z[num_atoms];
            for (i = 0; i < num_atoms; i++)
            {
                force_x[i] = 0;
                force_y[i] = 0;
                force_z[i] = 0;
            }
            get_energy(job, atomic, p, mol_nums, mol_lengths, bonds, num_bonds, num_atoms, force_x, force_y, force_z, &total_energy, &total_force);
            if (old_energy < total_energy)
            {
                hn = 0.2 * hn;
                old_force = total_force + 100;
            }
            else
            {
                hn = 1.2 * hn;
            }
            move(job, atomic, p, mol_nums, mol_lengths, num_atoms, force_x, force_y, force_z, hn);
            count_bad += 1;
            if (count_bad == 9999)
            {
                printf("danger");
            }
            //printf("%f %f %d\n", total_energy, total_force, p);
        }
        get_differences_restart(job, mol_nums, mol_lengths, differences, p, atomic);
    }
    update_pop(job, population, mol_nums, mol_lengths, atomic);
}

void update_gbest(struct Job job, struct PosAng gbest[job.molecule_types_number][job.max_mol_num],
                  struct PosAng population[job.population_per_core][job.molecule_types_number][job.max_mol_num],
                  int mol_lengths[job.molecule_types_number], double *gbest_chisq, int n_procs,
                  float all_chi_sq[n_procs][job.population_per_core], int rank, MPI_Comm comm,
                  struct Atom differences[job.population_per_core][job.molecule_types_number][job.max_mol_num][job.max_mol_length],
                  struct CharPair mol_nums[job.molecule_types_number], int iter, struct Bond bonds[job.max_mol_length], int num_bonds, struct Data exp_data[job.scattering_data_number][job.max_data_length], struct Data sim_data[job.scattering_data_number][job.max_data_length], int data_lengths[job.scattering_data_number])
{
    int best[2];
    float storage;
    float storage2;
    if (rank == 0){
        best_chi_sq(job, n_procs, all_chi_sq, best, gbest_chisq);
        printf("%d %d %d %f %f drew\n", iter, best[0], best[1], all_chi_sq[best[0]][best[1]], *gbest_chisq);
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
    if (storage2 <= storage && iter != 0)
    {
        if (rank == 0)
        {
            write_to_xyz(job, population, differences, mol_nums, mol_lengths, best[1], 0, iter);
            write_to_fit(job, population, differences, mol_nums, mol_lengths, best[1], 0, iter);
        }
    }
}

void get_atom_positions_restart(struct Job job,
                                struct Atom atomic[job.molecule_types_number][job.max_mol_num][job.max_mol_length],
                                struct CharPair mol_nums[job.molecule_types_number],
                                int mol_lengths[job.molecule_types_number])
{
    job.restart_file = fopen("restart.fit", "r");
    check_file_exists(job.restart_file, "restart.fit");
    char line[512];
    int line_count = 0;
    int mol_type_count = 0;
    int mol_num_count = 0;
    int mol_len_count = 0;
    while(fgets(line, sizeof(line), job.restart_file))
    {
        if(line_count > 1)
        {
            int j = 0;
            char *str;
            str = strtok(line, " ");
            while (str != NULL)
            {
                if (j == 1)
                {
                    strncpy(atomic[mol_type_count][mol_num_count][mol_len_count].label, str, 3);
                }
                if (j == 2)
                {
                    atomic[mol_type_count][mol_num_count][mol_len_count].x = atof(str);
                }
                if (j == 3)
                {
                    atomic[mol_type_count][mol_num_count][mol_len_count].y = atof(str);
                }
                if (j == 4)
                {
                    atomic[mol_type_count][mol_num_count][mol_len_count].z = atof(str);
                }
                int r;
                for (r = 0; r < job.scattering_data_number; r++)
                {
                    if (j == r + 5)
                    {
                        atomic[mol_type_count][mol_num_count][mol_len_count].sl[r] = atof(str);
                    }
                }
                j++;
                str = strtok(NULL, " ");
            }
            mol_len_count += 1;
            if (mol_len_count > mol_lengths[mol_type_count] - 1)
            {
                mol_len_count = 0;
                mol_num_count += 1;
            }
            if (mol_num_count > atoi(mol_nums[mol_type_count].keyword) - 1)
            {
                mol_len_count = 0;
                mol_num_count = 0;
                mol_type_count += 1;
            }
        }
        line_count += 1;
    }
    fclose(job.restart_file);
}



void get_bonds(struct Job job, struct CharPair mol_types[job.molecule_types_number], struct Bond bonds[job.max_mol_length], int *num_bonds)
{
    FILE *input_file;
    input_file = fopen("bonds.fit", "r");
    check_file_exists(input_file, "bonds.fit");
    char line[100];
    int i = 0;
    while(fgets(line, sizeof(line), input_file))
    {
        int j = 0;
        char *str;
        str = strtok(line, ";");
        while (str != NULL)
        {
            if (j == 0)
            {
                bonds[i].molecule_index = atoi(str);
            }
            if (j == 1)
            {
                bonds[i].a = atoi(str);
            }
            if (j == 2)
            {
                bonds[i].b = atoi(str);
            }
            j++;
            str = strtok(NULL, ";");
        }
        i++;
    }
    *num_bonds = i;
    fclose(input_file);
}


int main(int argc, char *argv[])
{
    int rank, n_procs;
    MPI_Comm comm = mpstart(&n_procs, &rank);
    double t1, t2;
    t1 = MPI_Wtime();
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
    struct Atom differences[job.population_per_core][job.molecule_types_number][job.max_mol_num][job.max_mol_length];
    struct PosAng population[job.population_per_core][job.molecule_types_number][job.max_mol_num];
    struct Bond bonds[job.max_mol_length];
    int num_bonds;

    get_atom_positions(job, molecules, mol_types);

    get_differences(job, mol_nums, mol_lengths, differences, molecules);

    get_bonds(job, mol_types, bonds, &num_bonds);
    /*else
    {
        struct Atom atomic[job.][job.molecule_types_number][job.max_mol_num][job.max_mol_length];

        get_atom_positions_restart(job, atomic, mol_nums, mol_lengths);

        get_differences_restart(job, mol_nums, mol_lengths, differences, 1, atomic);

        get_bonds(job, mol_types, bonds, &num_bonds);
    }*/
    build_population(job, population, mol_nums);

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
        energy_minimisation(job, population, differences, mol_nums, mol_lengths, bonds, num_bonds, exp_data, sim_data, data_lengths);
        float chi_sq[job.population_per_core];
        analyse(job, population, mol_nums, mol_lengths, differences, exp_data, sim_data, data_lengths, chi_sq);

        float all_chi_sq[n_procs][job.population_per_core];
        MPI_Gather(chi_sq, job.population_per_core, MPI_FLOAT, all_chi_sq, job.population_per_core, MPI_FLOAT, 0, comm);

        update_pbest(job, pbest, population, mol_lengths, pbest_chisq, chi_sq);
        update_gbest(job, gbest, population, mol_lengths, &gbest_chisq, n_procs, all_chi_sq, rank, comm, differences,
                     mol_nums, i, bonds, num_bonds, exp_data, sim_data, data_lengths);
        integrator(job, n_procs, rank, population, pbest, gbest, mol_lengths, velocity);
    }

    t2 = MPI_Wtime();
    if (rank == 0)
    {
        printf( "Elapsed time is %f\n", t2 - t1 );
    }
    MPI_Finalize();
}
