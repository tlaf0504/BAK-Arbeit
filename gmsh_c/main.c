#include <stdio.h>
#include "gmshc.h"

int main(int argc, char** argv) {
    int ierr;
    int argc_gmsh = 0;
    char* argv_gmsh[1] = {"/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/gmsh-4.3.0-Linux64/bin/gmsh"};

    gmshInitialize(argc_gmsh, argv_gmsh, 1, &ierr);
    gmshOptionSetNumber("General.Terminal", 1, &ierr);


    gmshOpen("/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/repo/problems/cylinder_cap/cylinder_cap.geo", &ierr);

    gmshModelMeshGenerate(2, &ierr);
    gmshModelMeshSetOrder(2, &ierr);

    gmshWrite("/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/repo/problems/cylinder_cap/cylinder_cap_api.msh22", &ierr);

    gmshFinalize(&ierr);






    return 0;
}