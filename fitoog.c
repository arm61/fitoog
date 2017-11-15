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
    float dq;
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
    strncpy(input[0].label, "name", sizeof("name"));
    strncpy(input[0].keyword, "fitoog_run", sizeof("fitoog_run"));
    strncpy(input[1].label, "num_molecules", sizeof("num_molecules"));
    strncpy(input[1].keyword, "1", sizeof("1"));
    strncpy(input[2].label, "restart", sizeof("restart"));
    strncpy(input[2].keyword, "no", sizeof("no"));
    strncpy(input[3].label, "cell", sizeof("name"));
    strncpy(input[3].keyword, "100", sizeof("100"));
    strncpy(input[4].label, "num_data", sizeof("num_data"));
    strncpy(input[4].keyword, "1", sizeof("1"));
    strncpy(input[5].label, "golden_vectors", sizeof("golden_vectors"));
    strncpy(input[5].keyword, "30", sizeof("30"));
    strncpy(input[6].label, "pop_size", sizeof("golden_vectors"));
    strncpy(input[6].keyword, "360", sizeof("30"));
    strncpy(input[7].label, "num_steps", sizeof("golden_vectors"));
    strncpy(input[7].keyword, "10", sizeof("30"));
}

void CheckFileExists(FILE *file, char *fileName)
{
    if (file == NULL)
    {
        printf("The %s file could not be found, as a result fitoog has quit.\n", fileName);
        exit(EXIT_FAILURE);
    }
}

char *FindWordInLine(char line[100], char word[50])
{
    char *found;
    found = strstr(line, word);
    return found;
}

void UpdateOutput(int numcPair, struct CharPair cPair[numcPair], char line[100], int i)
{
    const char wipe[50] = "                                                  ";
    strncpy(cPair[i].keyword, wipe, 50);
    strncpy(cPair[i].keyword, FirstCharAfterSpace(line), 50);
    strcat(cPair[i].keyword, "\0");
}

void FindUpdateOutput(int numcPair, struct CharPair cPair[numcPair], char line[100])
{
    int i;
    for (i = 0; i < numcPair; i++)
    {
        char *found = FindWordInLine(line, cPair[i].label);
        if (found != NULL)
        {
            UpdateOutput(numcPair, cPair, line, i);
        }
    }
}

void ReadInputFile(int numcPair, struct CharPair cPair[numcPair])
{
    FILE *inputfile;
    inputfile = fopen("fitoog.inp", "r");
    CheckFileExists(inputfile, "fitoog.inp");
    char line[100];
    while(fgets(line, sizeof(line), inputfile))
    {
        FindUpdateOutput(numcPair, cPair, line);
    }
    fclose(inputfile);
}

void MakeCell(int numcPair, struct CharPair cPair[numcPair], float cell[3])
{
    int i;
    for (i = 0; i < 3; i++)
    {
        cell[i] = atof(cPair[3].keyword);
    }
}

void CopyMoleculeLabels(char label[50], char id[50], int i)
{
    strncpy(label, id, 50);
    char str[50];
    sprintf(str, "%d", i);
    strcat(label, str);
    strcat(label, "\0");
}

void InitialiseMoleculeNames(int numMolecules, struct CharPair molTypeLoc[numMolecules], char id[50])
{
    int i;
    for (i = 0; i < numMolecules; i++)
    {
        CopyMoleculeLabels(molTypeLoc[i].label, id, i);
    }
}

void InitialiseAndRead(int numMolecules, struct CharPair molTypeLoc[numMolecules], char id[50])
{
    InitialiseMoleculeNames(numMolecules, molTypeLoc, id);
    ReadInputFile(numMolecules, molTypeLoc);
}

void GetFileLengths(int numFiles, int itemLengths[numFiles], struct CharPair itemType[numFiles])
{
    int i;
    for (i = 0; i < numFiles; i++)
    {
        FILE *inputFile;
        inputFile = fopen(itemType[i].keyword, "r");
        CheckFileExists(inputFile, itemType[i].keyword);
        char line[100];
        int j = 0;
        while (fgets(line, sizeof(line), inputFile))
        {
            j++;
        }
        itemLengths[i] = j;
        fclose(inputFile);
    }
}

int FindMaxInt(int length, int array[length])
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

int FindMaxCharPair(int length, struct CharPair array[length])
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

void ReadInDataInfo(char line[512], int numData, int maxDataLength, struct ExpData data[numData][maxDataLength], int k,
                    int i)
{
    int j = 0;
    int maxj = 0;
    char *str;
    str = strtok(line, ",");
    while (str != NULL)
    {
        if (j == 0)
        {
            data[i][k].q = atof(str);
        }
        if (j == 1)
        {
            data[i][k].i = atof(str);
        }
        if (j == 2)
        {
            data[i][k].di = atof(str);
        }
        if (j == 3)
        {
            data[i][k].dq = atof(str);
        }
        maxj = j;
        j++;
        str = strtok(NULL, ",");
    }
    if (maxj < 2)
    {
        data[i][k].di = data[i][k].i * 0.05;
        data[i][k].dq = data[i][k].q * 0.05;
    }
    if (maxj < 3)
    {
        data[i][k].dq = data[i][k].q * 0.05;
    }
}

void GetDataPoints(int numData, int maxDataLength, struct ExpData data[numData][maxDataLength],
                   struct CharPair dataTypes[numData])
{
    int i;
    for (i = 0; i < numData; i++)
    {
        FILE *dataFile;
        dataFile = fopen(dataTypes[i].keyword, "r");
        CheckFileExists(dataFile, dataTypes[i].keyword);
        char line[512];
        int k = 0;
        while (fgets(line, sizeof(line), dataFile))
        {
            ReadInDataInfo(line, numData, maxDataLength, data, k, i);
            k++;
        }
        fclose(dataFile);
    }
}

void ReadInAtomInfo(char line[512], int numMolecules, int maxMolLength, struct Atom readMols[numMolecules][maxMolLength], int k, int i, int numData)
{
    int j = 0;
    char *str;
    str strtok(line, " ");
    while (str != NULL)
    {
        if (j == 0)
        {
            strncpy(readMols[i][k].label, str, 3);
        }
        if (j == 1)
        {
            readMols[i][k].xpos = atof(str);
        }
        if (j == 2)
        {
            readMols[i][k].ypos = atof(str);
        }
        if (j == 3)
        {
            readMols[i][k].zpos = atof(str);
        }
        int r;
        for (r = 0; r < numData; r++)
        {
            if (j == r + 4)
            {
                readMols[i][k].scatLen[r] = atof(str);
            }
        }
        j++;
        str = strtok(NULL, " ");
    }
}

void GetAtomPositions(int numMolecules, int maxMolLength, struct Atom readMols[numMolecules][maxMolLength], struct CharPair molTypes[numMolecules], int numData)
{
    int i;
    for (i = 0; i < numMolecules; i++)
    {
        FILE *inputFile;
        inputFile = fopen(molTypes[i].keyword, "r");
        CheckFileExists(inputFile, molTypes[i].keyword);
        char line[512];
        int k = 0;
        while (fgets(line, sizeof(line), inputFile))
        {
            ReadInAtomInfo(line, numMolecules, maxMolLength, readMols, k, i, numData);
            k++;
        }
        fclose(inputFile);
    }
}

MPI_Comm mpstart(int *nProcs, int *rank)
{
    MPI_Init(NULL, NULL);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, nProcs);
    MPI_Comm_rank(comm, rank);
    return comm;
}

int main(int argc, char *argv[])
{
    int rank, nProcs;
    MPI_Comm comm = mpstart(&nProcs, &rank);
    srand(rank + 1 * time(NULL));

    int numInput = 7;
    struct CharPair input[numInput];

    initialiseKeywords(numInput, input);
    ReadInputFile(numInput, input);

    struct Job jobDetails;
    InitialiseJob(&jobDetails, numInput, input, nProcs);

    float cell[3];
    MakeCell(numInput, input, cell);

    int goldenVectors = jobDetails.goldenVectors;

    int numMolecules = jobDetails.numberMoleculeTypes;
    struct CharPair molTypes[numMolecules];
    InitialiseAndRead(numMolecules, molTypes, "molecule");

    int molLengths[numMolecules];
    GetFileLengths(numMolecules, molLengths, molTypes);
    int maxMolLength = FindMaxInt(numMolecules, molLengths);

    struct CharPair molNums[numMolecules];
    InitialiseAndRead(numMolecules, molNums, "num_mole");
    int maxMolNum = FindMaxCharPair(numMolecules, molNums);

    int numData = jobDetails.numberScatteringData;
    struct CharPair dataTypes[numData];
    InitialiseAndRead(numData, dataTypes, "data");

    int dataLengths[numData];
    GetFileLengths(numData, dataLengths, dataTypes);
    int maxDataLength = FindMaxInt(numData, dataLengths);

    struct ExpData data[numData][maxDataLength];
    struct ExpData simData[numData][maxDataLength];
    GetDataPoints(numData, maxDataLength, data, dataTypes);

    struct Atom readMols[numMolecules][maxMolLength];
    GetAtomPositions(numMolecules, maxMolLength, readMols, molTypes, numData);

    printf("Successful\n");
    MPI_Finalize();
}