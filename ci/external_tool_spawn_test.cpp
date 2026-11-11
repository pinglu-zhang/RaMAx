#include "external_tool.h"
#include "spawn_compat.h"
#include <cassert>
#include <csignal>
#include <fstream>
#include <iostream>
#include <thread>
#include <utility>
#include <sys/prctl.h>
#include <sys/wait.h>

namespace {
bool group = false;
int fail_init = 0, fail_flags = 0, fail_action = 0, live = 0;
size_t fdCount() {
    return std::distance(std::filesystem::directory_iterator("/proc/self/fd"),
                         std::filesystem::directory_iterator{});
}
}
extern "C" {
int __real_posix_spawn(pid_t*, const char*, const posix_spawn_file_actions_t*,
                      const posix_spawnattr_t*, char* const[], char* const[]);
int __real_posix_spawnattr_init(posix_spawnattr_t*);
int __real_posix_spawnattr_setflags(posix_spawnattr_t*, short);
int __real_posix_spawnattr_destroy(posix_spawnattr_t*);
int __real_posix_spawn_file_actions_adddup2(posix_spawn_file_actions_t*, int, int);
int __wrap_posix_spawnattr_init(posix_spawnattr_t* a) {
    if (std::exchange(fail_init, 0)) return ENOMEM;
    int r = __real_posix_spawnattr_init(a); if (!r) ++live; return r;
}
int __wrap_posix_spawnattr_setflags(posix_spawnattr_t* a, short f) {
    if (std::exchange(fail_flags, 0)) return EINVAL;
    return __real_posix_spawnattr_setflags(a, f);
}
int __wrap_posix_spawnattr_destroy(posix_spawnattr_t* a) {
    --live; return __real_posix_spawnattr_destroy(a);
}
int __wrap_posix_spawn_file_actions_adddup2(posix_spawn_file_actions_t* a, int f, int t) {
    if (std::exchange(fail_action, 0)) return EBADF;
    return __real_posix_spawn_file_actions_adddup2(a, f, t);
}
int __wrap_posix_spawn(pid_t* p, const char* e, const posix_spawn_file_actions_t* a,
                      const posix_spawnattr_t* s, char* const v[], char* const n[]) {
    short expected = group ? POSIX_SPAWN_SETPGROUP : 0;
#ifdef RAMAX_LEGACY_SPAWN
    expected |= POSIX_SPAWN_USEVFORK;
#endif
    short actual = 0;
    if (expected) {
        assert(s && posix_spawnattr_getflags(s, &actual) == 0 && actual == expected);
        if (group) { pid_t pg; assert(posix_spawnattr_getpgroup(s, &pg) == 0 && pg == 0); }
    } else assert(!s);
    return __real_posix_spawn(p, e, a, s, v, n);
}
}
int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "tree") {
        assert(getpgrp() == getpid());
        signal(SIGTERM, SIG_IGN);
        if (fork() == 0) { for (;;) pause(); }
        for (;;) pause();
    }
    assert(argc == 4);
    assert(prctl(PR_SET_CHILD_SUBREAPER, 1) == 0);
    std::filesystem::path root(argv[1]);
    std::filesystem::create_directories(root);
    auto out = root / "out", err = root / "err";
    const size_t fds = fdCount();
    auto r = RaMAxExternalTool::run("/bin/sh", {"-c", "printf stdout; printf stderr >&2; exit 7"}, out, err);
    assert(r.exit_code == 7 && !r.timed_out && !r.cancelled);
    assert(RaMAxExternalTool::readText(out) == "stdout");
    assert(RaMAxExternalTool::readText(err) == "stderr");
    r = RaMAxExternalTool::run("/bin/sh", {"-c", "kill -TERM $$"}, out, err);
    assert(r.termination_signal == SIGTERM && r.exit_code != 0);
    // Shared by Mash (no group) and wfmash/samtools (with group).
    for (bool grouped : {false, true}) {
        group = grouped;
        RaMAxExternalTool::RunOptions options; options.create_process_group = group;
        for (int* failure : {&fail_init, &fail_flags, &fail_action}) {
#ifndef RAMAX_LEGACY_SPAWN
            if (!group && failure != &fail_action) continue;
#endif
            *failure = 1;
            bool threw = false;
            try { RaMAxExternalTool::run("/bin/true", {}, out, err, options); }
            catch (const std::system_error&) { threw = true; }
            assert(threw && live == 0 && fdCount() == fds);
        }
    }
    group = true;
    for (bool cancel : {false, true}) {
        RaMAxExternalTool::RunOptions options;
        options.create_process_group = true;
        options.timeout = std::chrono::milliseconds(150);
        options.termination_grace = std::chrono::milliseconds(50);
        options.poll_interval = std::chrono::milliseconds(10);
        std::atomic<bool> requested{false};
        std::thread trigger;
        if (cancel) {
            options.timeout = std::chrono::milliseconds(0);
            options.cancellation_requested = &requested;
            trigger = std::thread([&] { std::this_thread::sleep_for(std::chrono::milliseconds(150)); requested = true; });
        }
        r = RaMAxExternalTool::run(std::filesystem::canonical("/proc/self/exe"), {"tree"}, out, err, options);
        if (trigger.joinable()) trigger.join();
        assert(r.cancelled == cancel && r.timed_out == !cancel && r.termination_signal == SIGKILL);
        int status;
        pid_t child = waitpid(-1, &status, WNOHANG);
        assert(child > 0 && WIFSIGNALED(status) && WTERMSIG(status) == SIGKILL);
        assert(waitpid(-1, &status, WNOHANG) == -1 && errno == ECHILD);
    }
    assert(live == 0 && fdCount() == fds);
    group = false;
    bool missing_threw = false;
    try {
        const auto missing = RaMAxExternalTool::run("/missing/executable", {}, out, err);
        // POSIX permits spawn to succeed and report exec failure as child exit
        // 127. glibc 2.17's vfork path uses this instead of returning ENOENT.
        missing_threw = missing.exit_code == 127;
    }
    catch (const std::runtime_error&) { missing_threw = true; }
    assert(missing_threw && live == 0 && fdCount() == fds);
    pid_t child = fork();
    assert(child >= 0);
    if (!child) {
        close(STDOUT_FILENO); close(STDERR_FILENO);
        auto closed = RaMAxExternalTool::run("/bin/sh",
            {"-c", "printf out; printf err >&2"}, out, err);
        assert(closed.exit_code == 0 && RaMAxExternalTool::readText(out) == "out" &&
               RaMAxExternalTool::readText(err) == "err");
        _exit(0);
    }
    int status;
    assert(waitpid(child, &status, 0) == child && WIFEXITED(status) && WEXITSTATUS(status) == 0);
    std::cout << "PASS: redirection, exits, signals, combined flags, failure cleanup, timeout/cancel descendants\n";
}
