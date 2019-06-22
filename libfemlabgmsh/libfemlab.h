//
// Created by tobiaslafer on 21.06.19.
//

#ifndef LIBFEMLABGMSH_LIBFEMLAB_H
#define LIBFEMLABGMSH_LIBFEMLAB_H

#include <stdlib.h>

void femlabGmshInitialize(int* ierr);
void femlabGmshFinalize(int* ierr);
void femlabGmshOpen(const char* filename, int* ierr);
void femlabGmshModelMeshSetOrder(const int order, int* ierr);
void femlabGmshModelMeshGenerate(const int dim, int* ierr);
void femlabGmshEnableConsoleLog(int* ierr);
void femlabGmshLoggerStart(int* ierr);
void femlabGmshLoggerStop(int* ierr);
void femlabGmshLoggerGet(char** log, int* ierr);




#endif //LIBFEMLABGMSH_LIBFEMLAB_H
