/*
 * Old libc has no posix_spawn chdir action. Spawn this small executable first,
 * then chdir and exec in its own address space, never in a vfork callback.
 * exec preserves PID, descriptors and process group for cleanup/cancellation.
 */
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: ramax-spawn-helper DIRECTORY EXECUTABLE [ARG ...]\n");
        return 126;
    }
    if (chdir(argv[1]) != 0) {
        fprintf(stderr, "ramax-spawn-helper: chdir %s: %s\n", argv[1], strerror(errno));
        return 126;
    }
    execvp(argv[2], argv + 2);
    int error = errno;
    fprintf(stderr, "ramax-spawn-helper: exec %s: %s\n", argv[2], strerror(error));
    return error == ENOENT ? 127 : 126;
}
