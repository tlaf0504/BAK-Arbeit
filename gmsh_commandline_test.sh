#!/bin/bash

GMSH_PATH=/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/gmsh-4.3.0-Linux64/bin
PROBLEM_PATH=/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/repo/problems


$GMSH_PATH/gmsh $PROBLEM_PATH/cylinder_cap/cylinder_cap.geo -2 -order 2 -format msh22
#$GMSH_PATH/gmsh $PROBLEM_PATH/cylinder_cap/cylinder_cap.geo -2
