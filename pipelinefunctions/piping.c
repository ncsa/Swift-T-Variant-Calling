/* pipe function definition */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void runpipe(int pfd[], char **cmd1, char **cmd2) {
     	int pid;
	switch (pid = fork()) {
	       	case 0: /* child */
			dup2(pfd[0], 0);
			close(pfd[1]); /* the child does not need this end of the pipe */
			execvp(*cmd2, cmd2);
			perror(*cmd2);
		default: /* parent */
			dup2(pfd[1], 1);
			close(pfd[0]); /* the parent does not need this end of the pipe */
			execvp(*cmd1, cmd1);
			perror(*cmd1);
		case -1:
			perror("fork");
			exit(1);
	}
}

