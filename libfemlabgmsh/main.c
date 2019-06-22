#include <stdio.h>
#include "libfemlab.h"

int main()
{
    int ierr;
    femlabGmshInitialize(&ierr);
    femlabGmshEnableConsoleLog(&ierr);
    femlabGmshLoggerStart(&ierr);


    femlabGmshOpen("/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/repo/problems/cylinder_cap/cylinder_cap.geo", &ierr);

    femlabGmshModelMeshGenerate(2, &ierr);
    femlabGmshModelMeshSetOrder(2, &ierr);



    char* log;
    femlabGmshLoggerGet(&log, &ierr);
    printf("%s", log);

    femlabGmshLoggerStop(&ierr);
    femlabGmshFinalize(&ierr);
}