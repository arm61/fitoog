//
// Created by arm61 on 14/11/17.
//

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <math.h>

#define PI 3.141592653
#define PHI 1.619033989

struct Job
{
    char name[50];
    int numberMoleculeTypes;
    int restart;
    float cell;
    int numberScatteringData;
    int goldenVectors;
    int populationSize;
    int numberSteps;
    int populationPerCore;
};

struct CharPair
{
    char label[50];
    char keyword[50];
};

struct Atom
{
    char label[3];
    float xpos;
    float ypos;
    float zpos;
    float scatLen[10];
};

struct PosAng
{
    float position[3];
    float angle[3];
};

struct ExpData
{
    float q;
    float i;
    float di;
};

char *FirstCharAfterSpace(char *input)
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

int TrueFalse(char *str)
{
    if (str == "yes")
    {
        return 1;
    }
    if (str == "no")
    {
        return 0;
    }
}

int CheckPopSizeVsNProcs(int populationSize, int nProcs)
{
    if (populationSize % nProcs != 0)
    {
        printf("For peak efficiency please ensure that the population can be spread evenly across the cores.\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    else {
        return populationSize / nProcs;
    }
}

void InitialiseJob(struct Job *jobDetails, int numInputs, struct CharPair input[numInputs], int nProcs)
{
    strncpy(jobDetails->name, input[0].keyword, 50);
    jobDetails->numberMoleculeTypes = atoi(input[1].keyword);
    jobDetails->restart = TrueFalse(input[2].keyword);
    jobDetails->cell = atof(input[3].keyword);
    jobDetails->numberScatteringData = atoi(input[4].keyword);
    jobDetails->goldenVectors = atoi(input[5].keyword);
    jobDetails->populationSize = atoi(input[6].keyword);
    jobDetails->numberSteps = atoi(input[7].keyword);
    jobDetails->populationPerCore = CheckPopSizeVsNProcs(jobDetails->populationSize, nProcs);
}

void initialiseKeywords(int numInput, struct CharPair input[numInput])
{
    strncpy(input[0].Label, "name", sizeof("name"));
    strncpy(input[0].Keyword, "fitoog_run", sizeof("fitoog_run"));
    strncpy(input[1].Label, "num_molecules", sizeof("num_molecules"));
    strncpy(input[1].Keyword, "1", sizeof("1"));
    strncpy(input[2].Label, "restart", sizeof("restart"));
    strncpy(input[2].Keyword, "no", sizeof("no"));
    strncpy(input[3].Label, "cell", sizeof("name"));
    strncpy(input[3].Keyword, "100", sizeof("100"));
    strncpy(input[4].Label, "num_data", sizeof("num_data"));
    strncpy(input[4].Keyword, "1", sizeof("1"));
    strncpy(input[5].Label, "golden_vectors", sizeof("golden_vectors"));
    strncpy(input[5].Keyword, "30", sizeof("30"));
    strncpy(input[6].Label, "pop_size", sizeof("golden_vectors"));
    strncpy(input[6].Keyword, "360", sizeof("30"));
    strncpy(input[7].Label, "num_steps", sizeof("golden_vectors"));
    strncpy(input[7].Keyword, "10", sizeof("30"));
}

int main(int argc, char *argv[])
{
    int rank, nProcs;
    MPI_Comm comm = mpstart(&nProcs, &rank);
    srand(rank + 1 * time(NULL));

    int numInput = 7;
    struct CharPair input[numInput];

    InitialiseKeywords(numInput, input);
}