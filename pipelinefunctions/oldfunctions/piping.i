%module pipingmodule
%{
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
%}

void runpipe(int pfd[], char **cmd1, char **cmd2);
