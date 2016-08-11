/* pipe demo */ /* Paul Krzyzanowski */
// from https://www.cs.rutgers.edu/~pxk/416/notes/c-tutorials/pipe.html
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define OK       0
#define NO_INPUT 1
#define TOO_LONG 2
#define CMD_LENGTH 30

void runpipe();
static int getLine (char *buff, size_t sz) ;

int main(int argc, char **argv)
{
	char ** cmd1;
	char ** cmd2;
	if(argc == 10){
		cmd1 = { argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8] };
		cmd2 = &argv[9];
	}
	int pid, status;
	int fd[2];
// reading user's commands:



	pipe(fd);
	switch (pid = fork()) {
		case 0: /* child */
			runpipe(fd, cmd1, cmd2);
			exit(0);
		default: /* parent */
			while ((pid = wait(&status)) != -1)
				fprintf(stderr, "process %d exits with %d\n", pid, WEXITSTATUS(status));
			break;
		case -1:
			perror("fork");
			exit(1);
	}
	exit(0);
}


void runpipe(int pfd[], char *cmd1[], char *cmd2[]) {
     	int pid;
	switch (pid = fork()) {
	       	case 0: /* child */
			dup2(pfd[0], 0);
			close(pfd[1]); /* the child does not need this end of the pipe */
			execvp(cmd2[0], cmd2);
			perror(cmd2[0]);
		default: /* parent */
			dup2(pfd[1], 1);
			close(pfd[0]); /* the parent does not need this end of the pipe */
			execvp(cmd1[0], cmd1);
			perror(cmd1[0]);
		case -1:
//			perror("fork");
			exit(1);
	}
}

static int getLine (char *buff, size_t sz) {
            int ch, extra;
            if (fgets (buff, sz, stdin) == NULL)
                return NO_INPUT;
// If it was too long, there'll be no newline. In that case, we flush
// to end of line so that excess doesn't affect the next call.
            if (buff[strlen(buff)-1] != '\n') {
                extra = 0;
                while (((ch = getchar()) != '\n') && (ch != EOF))
                        extra = 1;
                return (extra == 1) ? TOO_LONG : OK;
            }
// Otherwise remove newline and give string back to caller.
             buff[strlen(buff)-1] = '\0';
             return OK;
}
