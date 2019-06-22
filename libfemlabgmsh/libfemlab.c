//
// Created by tobiaslafer on 21.06.19.
//

#include "libfemlab.h"
#include "gmshc.h"
#include <string.h>


void femlabGmshInitialize(int* ierr)
{
    int argc = 0;
    gmshInitialize(argc, NULL, 1, ierr);
}

void femlabGmshFinalize(int* ierr)
{
    gmshFinalize(ierr);
}

void femlabGmshClear(int* ierr)
{
    gmshClear(ierr);
}

void femlabGmshWrite(const char* filename, int* ierr)
{
    gmshWrite(filename, ierr);
}

void femlabGmshOpen(const char* filename, int* ierr)
{
    gmshOpen(filename, ierr);
}

void femlabGmshMerge(const char* filename, int* ierr)
{
    gmshMerge(filename, ierr);
}

void femlabGmshModelMeshRefine(int* ierr)
{
    gmshModelMeshRefine(ierr);
}

void femlabGmshModelMeshSetOrder(const int order, int* ierr)
{
    gmshModelMeshSetOrder(order, ierr);
}


void femlabGmshModelMeshGenerate(const int dim, int* ierr)
{
    gmshModelMeshGenerate(dim, ierr);
}


void femlabGmshLoggerStart(int* ierr)
{
    gmshLoggerStart(ierr);
}

void femlabGmshLoggerStop(int* ierr)
{
    gmshLoggerStop(ierr);
}

void femlabGmshEnableConsoleLog(int* ierr)
{
    gmshOptionSetNumber("General.Terminal", 1, ierr);
}

void femlabGmshLoggerGet(char** log, int* ierr)
{
    char** log_array;
    size_t n_log_entries;

    gmshLoggerGet(&log_array, &n_log_entries, ierr);

    unsigned int log_string_length = 0;
    for (int i = 0; i < n_log_entries; i++)
    {
        char* current_entry = log_array[i];
        log_string_length += (strlen(current_entry) + 1);
    }

    // ATTENTION: Provide mechanism for freeing memory after testing the basic functionality
    *log = malloc(log_string_length  + 1);
    if (*log == NULL)
    {
        *ierr = -1;
        return;
    }

    memset(*log, 0, log_string_length + 1);

    for (int i = 0; i < n_log_entries; i++)
    {
        char* current_entry = log_array[i];
        unsigned int current_length = strlen(current_entry);

        // Temporary string containing log entry plus a '\n'
        char* tmp = malloc(current_length + 2);
        if (tmp == NULL)
        {
            *ierr = -1;
            return;
        }

        memset(tmp, 0, current_length + 2);
        strcpy(tmp, log_array[i]);
        strcat(tmp, "\n");

        strcat(*log, tmp);
        free(tmp);
    }
}

