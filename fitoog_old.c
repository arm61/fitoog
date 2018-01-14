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
    str = strtok(line, " ");
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
    str = strtok(line, " ");
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

void PopulateDifferences(int numMolecules, int maxMolNum, struct CharPair molNums[numMolecules],
                         int maxMolLength, int molLengths[numMolecules],
                         struct Atom differences[numMolecules][maxMolNum][maxMolLength],
                         struct Atom readMols[numMolecules][maxMolLength], int numData)
{
    float mean[numMolecules][3];
    int i;
    for (i = 0; i < numMolecules; i++)
    {
        int j;
        for (j = 0; j < 3; j++)
        {
            mean[i][j] = 0.;
        }
        int k;
        for (k = 0; k < molLengths[i]; k++)
        {
            mean[i][0] += readMols[i][k].xpos;
            mean[i][1] += readMols[i][k].ypos;
            mean[i][2] += readMols[i][k].zpos;
        }
        for (k = 0; k < 3; k++)
        {
            mean[i][k] /= molLengths[i];
        }
    }
    for (i = 0; i < numMolecules; i++)
    {
        int j;
        for (j = 0; j < atoi(molNums[i].keyword); j++)
        {
            int k;
            for (k = 0; k < molLengths[i]; k++)
            {
                strncpy(differences[i][j][k].label, readMols[i][k].label, 3);
                differences[i][j][k].xpos = readMols[i][k].xpos - mean[i][0];
                differences[i][j][k].ypos = readMols[i][k].ypos - mean[i][1];
                differences[i][j][k].zpos = readMols[i][k].zpos - mean[i][2];
                int r = 0;
                for (r = 0; r < numData; r++)
                {
                    differences[i][j][k].scatLen[r] = readMols[i][k].scatLen[r];
                }
            }
        }
    }
}

float RandFloat(float max, float min)
{
    float x = (float)rand() / ((float)RAND_MAX / max) + min;
    return x;
}

void MakeRandom(int numMolecules, int maxMolNum, struct CharPair molNums[numMolecules],
                struct PosAng positions[numMolecules][maxMolNum], float cell[3])
{
    int i;
    for (i = 0; i < numMolecules; i++)
    {
        int j;
        for (j = 0; j < atoi(molNums[i].keyword); j++)
        {
            int k;
            for (k = 0; k < 3; k++)
            {
                positions[i][j].position[k] = RandFloat(cell[k], 0);
                positions[i][j].angle[k] = RandFloat(2 * PI, 0);
            }
        }
    }
}

void MakePopulation(int populationSize, int numMolecules, int maxMolNum,
                    struct PosAng population[populationSize][numMolecules][maxMolNum],
                    struct CharPair molNums[numMolecules], float cell[3])
{
    int i;
    for (i = 0; i < populationSize; i++)
    {
        MakeRandom(numMolecules, maxMolNum, molNums, population[i], cell);
    }
}

void RotationMatrix(float x, float y, float z, float rm[3][3]){
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

void ConvertToAtomic(int numMolecules, int maxMolNum, struct CharPair molNums[numMolecules], int maxMolLength,
                     int molLengths[numMolecules], struct Atom differences[numMolecules][maxMolNum][maxMolLength],
                     struct Atom atomic[numMolecules][maxMolNum][maxMolLength],
                     struct PosAng positions[numMolecules][maxMolNum], int numData){
    int i;
    for (i = 0; i < numMolecules; i++){
        int j;
        for (j = 0; j < atoi(molNums[i].keyword); j++){
            float rm[3][3];
            RotationMatrix(positions[i][j].angle[0], positions[i][j].angle[1], positions[i][j].angle[2], rm);
            int k;
            for (k = 0; k < molLengths[i]; k++){
                strncpy(atomic[i][j][k].label, differences[i][j][k].label, 3);
                atomic[i][j][k].xpos = positions[i][j].position[0] - differences[i][j][k].xpos;
                atomic[i][j][k].ypos = positions[i][j].position[1] - differences[i][j][k].ypos;
                atomic[i][j][k].zpos = positions[i][j].position[2] - differences[i][j][k].zpos;
                int r = 0;
                for (r = 0; r < numData; r++){
                    atomic[i][j][k].scatLen[r] = differences[i][j][k].scatLen[r];
                }
            }
            struct Atom other[molLengths[i]];
            for (k = 0; k < molLengths[i]; k++) {
                other[k].xpos = atomic[i][j][k].xpos - positions[i][j].position[0];
                other[k].ypos = atomic[i][j][k].ypos - positions[i][j].position[1];
                other[k].zpos = atomic[i][j][k].zpos - positions[i][j].position[2];
            }
            struct Atom another[molLengths[i]];
            for (k = 0; k < molLengths[i]; k++){
                another[k].xpos = rm[0][0] * other[k].xpos + rm[0][1] * other[k].ypos + rm[0][2] * other[k].zpos;
                another[k].ypos = rm[1][0] * other[k].xpos + rm[1][1] * other[k].ypos + rm[1][2] * other[k].zpos;
                another[k].ypos = rm[2][0] * other[k].xpos + rm[2][1] * other[k].ypos + rm[2][2] * other[k].zpos;
            }
            for (k = 0; k < molLengths[i]; k++){
                atomic[i][j][k].xpos = another[k].xpos + positions[i][j].position[0];
                atomic[i][j][k].ypos = another[k].ypos + positions[i][j].position[1];
                atomic[i][j][k].zpos = another[k].zpos + positions[i][j].position[2];
            }
        }
    }
}

void Normalise(int numData, int maxDataLength, struct ExpData data[numData][maxDataLength],
               struct ExpData simData[numData][maxDataLength], int dataLengths[numData], int i, int rank)
{
    int random = rand() % dataLengths[i];
    float norm = data[i][random].i / simData[i][random].i;
    int j = 0;
    for (j = 0; j < dataLengths[i]; j++)
    {
        simData[i][j].i = simData[i][j].i * norm;
    }
}

void DiffractionCalculator(int numMolecules, int maxMolNum, int maxMolLength,
                           struct Atom atomic [numMolecules][maxMolNum][maxMolLength],
                           struct CharPair molNums[numMolecules], int molLengths[numMolecules], int numData,
                           int maxDataLength, struct ExpData data[numData][maxDataLength], int dataLengths[numData],
                           int goldenVectors, struct ExpData simData[numData][maxDataLength], int rank){
    int i = 0;
    for (i = 0; i < numData; i++){
        int q;
        for (q = 0; q < dataLengths[i]; q++){
            float sum_int = 0.;
            float k_min = (-1. * ((float)goldenVectors - 1.)) / 2.;
            float k_max = ((float)goldenVectors - 1.) / 2.;
            float a;
            for (a = k_min; a < (k_max + 1); a += 1){
                float qx = data[i][q].q * cos(asin((2. * a) / (float)goldenVectors)) * cos((2. * PI * a) / (PHI));
                float qy = data[i][q].q * cos(asin((2. * a) / (float)goldenVectors)) * sin((2. * PI * a) / (PHI));
                float qz = (2.0 * a * data[i][q].q) / (float)goldenVectors;
                float sum_coss = 0.0;
                float sum_sinn = 0.0;
                int j;
                for (j = 0; j < numMolecules; j++){
                    int k;
                    for (k = 0; k < atoi(molNums[j].keyword); k++){
                        int l;
                        for (l = 0; l < molLengths[j]; l++){
                            float dot = qx * atomic[j][k][l].xpos +
                                    qy * atomic[j][k][l].ypos +
                                    qz * atomic[j][k][l].zpos;
                            float coss = atomic[j][k][l].scatLen[i] * cos(dot);
                            float sinn = atomic[j][k][l].scatLen[i] * sin(dot);
                            sum_coss += coss;
                            sum_sinn += sinn;
                        }
                        float sq_coss = sum_coss * sum_coss;
                        float sq_sinn = sum_sinn * sum_sinn;
                        sum_int += (sq_coss + sq_sinn);
                    }
                }
            }
            float inten = sum_int / goldenVectors;
            simData[i][q].i = inten;
            simData[i][q].q = data[i][q].q;
        }
        Normalise(numData, maxDataLength, data, simData, dataLengths, i, rank);
    }
}

void ChiSquared(int numData, int maxDataLength, struct ExpData data[numData][maxDataLength], int dataLengths[numData],
                struct ExpData simData[numData][maxDataLength], float chiSq[numData], float cell[3])
{
    int i;
    for (i = 0; i < numData; i++)
    {
        chiSq[i] = 0.;
        int j;
        for (j = 0; j < dataLengths[i]; j++)
        {
            if (data[i][j].q > (2 * PI) / cell[0])
            {
                chiSq[i] += ((simData[i][j].i - data[i][j].i) * (simData[i][j].i - data[i][j].i)) /
                        (data[i][j].di * data[i][j].di);
            }
        }
    }
}

void Analysis(int populationPerCore, int numMolecules, int maxMolNum, int maxMolLength,
              struct CharPair molNums[numMolecules], int molLengths[numMolecules],
              struct Atom differences[numMolecules][maxMolNum][maxMolLength],
              struct PosAng population[populationPerCore][numMolecules][maxMolNum], int numData, int maxDataLength,
              struct ExpData data[numData][maxDataLength], struct ExpData simData[numData][maxDataLength],
              int dataLengths[numData], int goldenVectors, int rank, float chiSq[populationPerCore][numData],
              float cell[3])
{
    int i;
    for (i  = 0; i < populationPerCore; i++)
    {
        struct Atom atomic[numMolecules][maxMolNum][maxMolLength];
        ConvertToAtomic(numMolecules, maxMolNum, molNums, maxMolLength, molLengths, differences, atomic, population[i],
                        numData);
        DiffractionCalculator(numMolecules, maxMolNum, maxMolLength, atomic, molNums, molLengths, numData,
                              maxDataLength, data, dataLengths, goldenVectors, simData, rank);
        ChiSquared(numData, maxDataLength, data, dataLengths, simData, chiSq[i], cell);
    }
}

void flatten3(int a, int b, int c, float array[a][b][c], float flat[a*b*c])
{
    int i;
    for (i = 0; i < a; i++)
    {
        int j;
        for (j = 0;  j < b; j++)
        {
            int k;
            for (k = 0; k < c; k++)
            {
                int comb = (i + a * (j + b * k));
                flat[comb] = array[i][j][k];
            }
        }
    }
}

void unflatten3(int a, int b, int c, float array[a][b][c], float flat[a*b*c])
{
    int comb = 0;
    int i;
    for (i = 0; i < a; i++)
    {
        int j;
        for (j = 0;  j < b; j++)
        {
            int k;
            for (k = 0; k < c; k++)
            {
                array[i][j][k] = flat[comb];
                comb++;
            }
        }
    }
}

void flatten3PosAng(int a, int b, int c, int d, int e, struct PosAng array[a][b][c][d][e],
                    struct PosAng flat[a*b*c][d][e], int z[d])
{
    int i;
    for (i = 0; i < a; i++)
    {
        int j;
        for (j = 0; j < b; j++)
        {
            int k;
            for (k = 0; k < c; k++)
            {
                int comb = (i + a * (j + b * k));
                int l;
                for (l = 0; l < d; l++)
                {
                    int m;
                    for (m = 0; m < z[l]; m++)
                    {
                        flat[comb][l][m] = array[i][j][k][l][m];
                    }
                }
            }
        }
    }
}

void unflatten3PosAng(int a, int b, int c, int d, int e, struct PosAng array[a][b][c][d][e],
                   struct PosAng flat[a*b*c][d][e], int z[d])
{
    int i;
    for (i = 0; i < a; i++)
    {
        int j;
        for (j = 0; j < b; j++)
        {
            int k;
            for (k = 0; k < c; k++)
            {
                int comb = (i + a * (j + b * k));
                int l;
                for (l = 0; l < d; l++)
                {
                    int m;
                    for (m = 0; m < z[l]; m++)
                    {
                        array[i][j][k][l][m] = flat[comb][l][m];
                    }
                }
            }
        }
    }
}

void popAssign(int a, int b, struct PosAng z[a][b], struct PosAng y[a][b])
{
    int i;
    for (i = 0; i < a; i++)
    {
        int j;
        for (j = 0; j < b; j++)
        {
            int k;
            for (k = 0; k < 3; k++)
            {
                y[i][j].position[k] = z[i][j].position[k];
                y[i][j].angle[k] = z[i][j].angle[k];
            }
        }
    }
}

void quickSort(int nProcs, int populationPerCore, int numData, float allChiSq[nProcs * populationPerCore * numData],
               int numMolecules, int maxMolNum,
               struct PosAng allPopulation[nProcs * populationPerCore * numData][numMolecules][maxMolNum], int low,
               int high)
{
    int pivot, j, temp, i;
    struct PosAng tempPop[numMolecules][maxMolNum];
    if (low < high)
    {
        pivot = low;
        i = low;
        j = high;
        while (i < j)
        {
            while ((allChiSq[i] <= allChiSq[pivot]) && (i < high))
            {
                i++;
            }
            while (allChiSq[j] > allChiSq[pivot])
            {
                j--;
            }
            if (i < j)
            {
                temp = allChiSq[i];
                popAssign(numMolecules, maxMolNum, allPopulation[i], tempPop);
                allChiSq[i] = allChiSq[j];
                popAssign(numMolecules, maxMolNum, allPopulation[j], allPopulation[i]);
                allChiSq[j] = temp;
                popAssign(numMolecules, maxMolNum, tempPop, allPopulation[j]);
            }
        }
        temp = allChiSq[pivot];
        popAssign(numMolecules, maxMolNum, allPopulation[pivot], tempPop);
        allChiSq[pivot] = allChiSq[j];
        popAssign(numMolecules, maxMolNum, allPopulation[j], allPopulation[pivot]);
        allChiSq[j] = temp;
        popAssign(numMolecules, maxMolNum, tempPop, allPopulation[j]);
        quickSort(nProcs, populationPerCore, numData, allChiSq, numMolecules, maxMolNum, allPopulation, low, j-1);
        quickSort(nProcs, populationPerCore, numData, allChiSq, numMolecules, maxMolNum, allPopulation, j+1, high);
    }
}

void sort(int nProcs, int populationPerCore, int numData, float allChiSq[nProcs][populationPerCore][numData],
          int numMolecules, int maxMolNum,
          struct PosAng allPopulation[nProcs][populationPerCore][numData][numMolecules][maxMolNum],
          int molLengths[numMolecules])
{
    float allChiSqFlat[nProcs * populationPerCore * numData];
    flatten3(nProcs, populationPerCore, numData, allChiSq, allChiSqFlat);
    struct PosAng allPopulationFlat[nProcs * populationPerCore * numData][numMolecules][maxMolNum];
    flatten3PosAng(nProcs, populationPerCore, numData, numMolecules, maxMolNum, allPopulation, allPopulationFlat, molLengths);
    quickSort(nProcs, populationPerCore, numData, allChiSqFlat, numMolecules, maxMolNum, allPopulationFlat, 0, (nProcs *
            populationPerCore * numData)-1);
    unflatten3(nProcs, populationPerCore, numData, allChiSq, allChiSqFlat);
    unflatten3PosAng(nProcs, populationPerCore, numData, numMolecules, maxMolNum, allPopulation, allPopulationFlat, molLengths);
}

void transferPrep(int populationPerCore, int numData, int numMolecules, int maxMolNum,
                  struct PosAng population[populationPerCore][numMolecules][maxMolNum],
                  float populationTrans[populationPerCore][numData][numMolecules][maxMolNum][6],
                  int molLengths[numMolecules])
{
    int i;
    for (i = 0; i < populationPerCore; i++)
    {
        int m;
        for (m = 0; m < numData; m++)
        {
            int j;
            for (j = 0; j < numMolecules; j++)
            {
                int k;
                for (k = 0; k < molLengths[j]; k++)
                {
                    populationTrans[i][m][j][k][0] = population[i][j][k].position[0];
                    populationTrans[i][m][j][k][1] = population[i][j][k].position[1];
                    populationTrans[i][m][j][k][2] = population[i][j][k].position[2];
                    populationTrans[i][m][j][k][3] = population[i][j][k].angle[0];
                    populationTrans[i][m][j][k][4] = population[i][j][k].angle[1];
                    populationTrans[i][m][j][k][5] = population[i][j][k].angle[2];
                }
            }
        }
    }
}

void transferReturn(int nProcs, int populationPerCore, int numData, int numMolecules, int maxMolNum,
                    struct PosAng population[nProcs][populationPerCore][numData][numMolecules][maxMolNum],
                    float populationTrans[nProcs][populationPerCore][numData][numMolecules][maxMolNum][6],
                    int molLengths[numMolecules])
{
    int a;
    for (a = 0; a < nProcs; a++)
    {
        int i;
        for (i = 0; i < populationPerCore; i++)
        {
            int m;
            for (m = 0; m < numData; m++)
            {
                int j;
                for (j = 0; j < numMolecules; j++)
                {
                    int k;
                    for (k = 0; k < molLengths[j]; k++)
                    {
                        population[a][i][m][j][k].position[0] = populationTrans[a][i][m][j][k][0];
                        population[a][i][m][j][k].position[1] = populationTrans[a][i][m][j][k][1];
                        population[a][i][m][j][k].position[2] = populationTrans[a][i][m][j][k][2];
                        population[a][i][m][j][k].angle[0] = populationTrans[a][i][m][j][k][3];
                        population[a][i][m][j][k].angle[1] = populationTrans[a][i][m][j][k][4];
                        population[a][i][m][j][k].angle[2] = populationTrans[a][i][m][j][k][5];
                    }
                }
            }
        }
    }
}

void transferPopulations(int populationPerCore, int numData, int numMolecules, int maxMolNum,
                         struct PosAng population[populationPerCore][numMolecules][maxMolNum],
                         int molLengths[numMolecules], int nProcs, MPI_Comm comm,
                         struct PosAng allPopulation[nProcs][populationPerCore][numData][numMolecules][maxMolNum])
{
    float populationTrans[populationPerCore][numData][numMolecules][maxMolNum][6];
    transferPrep(populationPerCore, numData, numMolecules, maxMolNum, population, populationTrans,
                 molLengths);

    float allPopulationTrans[nProcs][populationPerCore][numData][numMolecules][maxMolNum][6];
    MPI_Gather(populationTrans, numData * populationPerCore * numData * numMolecules * maxMolNum * 6,
               MPI_FLOAT,
               allPopulationTrans, numData * populationPerCore * numData * numMolecules * maxMolNum * 6,
               MPI_FLOAT, 0, comm);

    transferReturn(nProcs, populationPerCore, numData, numMolecules, maxMolNum, allPopulation,
                   allPopulationTrans, molLengths);

}

void PrintToXYZ(int numMolecules, int maxMolNum, struct CharPair molNums[numMolecules], int maxMolLength,
                int molLengths[numMolecules], struct Atom atomic[numMolecules][maxMolNum][maxMolLength], int rank){
    int i;
    int atoms = 0;
    for (i = 0; i < numMolecules; i++){
        int j;
        for (j = 0; j < atoi(molNums[i].keyword); j++) {
            atoms += molLengths[i];
        }
    }
    FILE *f = fopen("output.xyz", "w");
    CheckFileExists(f, "output.xyz");
    fprintf(f, "%d\n", atoms);
    fprintf(f, "produced by arm61\n");
    for (i = 0; i < numMolecules; i++){
        int j;
        for (j = 0; j < atoi(molNums[i].keyword); j++) {
            int k;
            for (k = 0; k < molLengths[i]; k++){
                fprintf(f, "%s %f %f %f\n", atomic[i][j][k].label, atomic[i][j][k].xpos, atomic[i][j][k].ypos,
                        atomic[i][j][k].zpos);
            }
        }
    }
    fclose(f);
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

    struct Atom differences[numMolecules][maxMolNum][maxMolLength];
    PopulateDifferences(numMolecules, maxMolNum, molNums, maxMolLength, molLengths, differences, readMols, numData);

    struct PosAng population[jobDetails.populationPerCore][numMolecules][maxMolNum];
    MakePopulation(jobDetails.populationPerCore, numMolecules, maxMolNum, population, molNums, cell);

    float chiSq[jobDetails.populationPerCore][numData];
    Analysis(jobDetails.populationPerCore, numMolecules, maxMolNum, maxMolLength, molNums, molLengths, differences,
             population, numData, maxDataLength, data, simData, dataLengths, goldenVectors, rank, chiSq, cell);

    float allChiSq[nProcs][jobDetails.populationPerCore][numData];
    MPI_Gather(chiSq, numData * jobDetails.populationPerCore, MPI_FLOAT, allChiSq,
               numData * jobDetails.populationPerCore, MPI_FLOAT, 0, comm);

    struct PosAng allPopulation[nProcs][jobDetails.populationPerCore][numData][numMolecules][maxMolNum];
    transferPopulations(jobDetails.populationPerCore, numData, numMolecules, maxMolNum, population, molLengths, nProcs,
                        comm, allPopulation);

    if (rank == 0){
        sort(nProcs, jobDetails.populationPerCore, numData, allChiSq, numMolecules, maxMolNum, allPopulation,
             molLengths);
        struct Atom bestAtomic[numMolecules][maxMolNum][maxMolLength];
        ConvertToAtomic(numMolecules, maxMolNum, molNums, maxMolLength, molLengths, differences, bestAtomic, allPopulation[0][0][0],
                        numData);
        PrintToXYZ(numMolecules, maxMolNum, molNums, maxMolLength, molLengths, bestAtomic, rank);
    }

    printf("Successful\n");
    MPI_Finalize();
}