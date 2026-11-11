#pragma once

#include <spawn.h>
#include <system_error>
#include <unistd.h>

// Test builds exercise the old-header branch without redefining libc macros.
#if defined(__linux__) && defined(__GLIBC__) && defined(POSIX_SPAWN_USEVFORK)
#if !__GLIBC_PREREQ(2, 24) || defined(RAMAX_TEST_LEGACY_SPAWN)
#define RAMAX_LEGACY_SPAWN 1
#endif
#endif
#if defined(__linux__) && defined(__GLIBC__)
#if __GLIBC_PREREQ(2, 29) && !defined(RAMAX_TEST_LEGACY_SPAWN)
#define RAMAX_SPAWN_HAS_CHDIR 1
#endif
#endif

namespace RaMAxSpawn {

inline void check(int error, const char* operation) {
    if (error) throw std::system_error(error, std::generic_category(), operation);
}

class Attributes {
public:
    explicit Attributes(bool process_group = false) {
        short flags = process_group ? POSIX_SPAWN_SETPGROUP : 0;
#ifdef RAMAX_LEGACY_SPAWN
        // libc owns the vfork child; no application code runs in the parent's
        // address space. Preserve process groups used by timeout/cancellation.
        flags |= POSIX_SPAWN_USEVFORK;
#endif
        if (!flags) return;
        check(posix_spawnattr_init(&value_), "Cannot initialize spawn attributes");
        initialized_ = true;
        int error = process_group ? posix_spawnattr_setpgroup(&value_, 0) : 0;
        if (!error) error = posix_spawnattr_setflags(&value_, flags);
        if (error) {
            posix_spawnattr_destroy(&value_);
            initialized_ = false;
            check(error, "Cannot configure spawn attributes");
        }
    }
    ~Attributes() { if (initialized_) posix_spawnattr_destroy(&value_); }
    Attributes(const Attributes&) = delete;
    Attributes& operator=(const Attributes&) = delete;
    const posix_spawnattr_t* get() const { return initialized_ ? &value_ : nullptr; }
private:
    posix_spawnattr_t value_{};
    bool initialized_ = false;
};

class FileActions {
public:
    FileActions() { check(posix_spawn_file_actions_init(&value_), "Cannot initialize spawn actions"); }
    ~FileActions() { posix_spawn_file_actions_destroy(&value_); }
    FileActions(const FileActions&) = delete;
    FileActions& operator=(const FileActions&) = delete;
    void dup2(int from, int to) {
        check(posix_spawn_file_actions_adddup2(&value_, from, to), "Cannot redirect child descriptor");
    }
    void close(int fd) {
        check(posix_spawn_file_actions_addclose(&value_, fd), "Cannot close child descriptor");
    }
#ifdef RAMAX_SPAWN_HAS_CHDIR
    void chdir(const char* path) {
        check(posix_spawn_file_actions_addchdir_np(&value_, path), "Cannot set child working directory");
    }
#endif
    const posix_spawn_file_actions_t* get() const { return &value_; }
private:
    posix_spawn_file_actions_t value_{};
};

class FileDescriptor {
public:
    explicit FileDescriptor(int fd) : fd_(fd) {}
    ~FileDescriptor() { if (fd_ >= 0) ::close(fd_); }
    FileDescriptor(const FileDescriptor&) = delete;
    FileDescriptor& operator=(const FileDescriptor&) = delete;
    int get() const { return fd_; }
private:
    int fd_;
};
} // namespace RaMAxSpawn
