/* pipe demo */ /* Paul Krzyzanowski */
// from https://www.cs.rutgers.edu/~pxk/416/notes/c-tutorials/pipe.html
//

// Call the executable as: 
// ./pipe <no. of arguments of the first command> <no. of arguments of the second command> <first command> <seperated list of arguments fort the first comannd> <second command> <seperated list of arguments for the second command>


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void runpipe(); 

int main(int argc, char *argv[] )
{
	int pid, status,i;
	int fd[2];

	// reading user's commands
	char *cmd1[ *(argv[1]) ];
	char *cmd2[ *(argv[2]) ];
	printf("\n number of parameters = %d", argc);
//	printf("\n\n cmd1 is %s \t\t", *cmd1 );
	printf("\nsize of cmd1 %s", argv[1] );
	int x =  *argv[1];
	printf("\t\t\t, x is %d\t\t\t\t\t\t\t\t\t\n",x);

//	printf("\n\n cmd2 is %s \t\t", *cmd2);
	printf("\nsize of cmd2 %s \n", argv[2] );

	for ( i = 3; i <= *argv[1]+3; i++ ) {
		cmd1[i-3] =  argv[i];
	}
//	printf("\n\n After initializaion, cmd1 is %s \t\t", *cmd1 );

	for ( i = *argv[1]+4; i < argc; i++ ) {
                 cmd2[i-(*argv[1])] =  argv[i];
	}
//	printf("\n\n After initializaion, cmd2 is %s \t\t", *cmd2);

	if ( *argv[1] < argc )  {
		exit(1);
	}
	
//	pipe(fd); 
//	switch (pid = fork()) { 
//		case 0: /* child */ 
//			runpipe(fd, &cmd1, &cmd2); 
//			exit(0); 
//		default: /* parent */ 
//			while ((pid = wait(&status)) != -1) 
//				fprintf(stderr, "process %d exits with %d\n", pid, WEXITSTATUS(status)); 
//			break; 
//		case -1: 
//			perror("fork"); 
//			exit(1); 
//	} 
//	exit(0);
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
			perror("fork");
			exit(1);
	}
}

