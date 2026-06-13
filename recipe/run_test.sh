#!/usr/bin/env bash
set -euxo pipefail

ramax_bin="${PREFIX}/bin/ramax"

if "${ramax_bin}" --help | grep -q -- "--mask-repeats"; then
  echo "--mask-repeats should not be exposed in packaged builds" >&2
  exit 1
fi

glibc_versions="$(
  objdump -T "${ramax_bin}" \
    | grep -o 'GLIBC_[0-9.]*' \
    | sort -Vu
)"
printf '%s\n' "${glibc_versions}"

max_glibc="$(printf '%s\n' "${glibc_versions}" | tail -n 1)"
test -n "${max_glibc}"
test "$(
  printf '%s\n' "${max_glibc}" GLIBC_2.19 | sort -V | tail -n 1
)" = "GLIBC_2.19"

ldd_output="$(ldd "${ramax_bin}")"
printf '%s\n' "${ldd_output}"

check_conda_library() {
  local library_pattern="$1"
  local line
  local resolved

  line="$(printf '%s\n' "${ldd_output}" | grep -E "[[:space:]]${library_pattern}[^[:space:]]*[[:space:]]+=>" | head -n 1)"
  test -n "${line}"

  resolved="$(printf '%s\n' "${line}" | awk '{print $3}')"
  case "${resolved}" in
    "${PREFIX}/"*) ;;
    *)
      echo "Expected ${library_pattern} to resolve inside ${PREFIX}, got ${resolved}" >&2
      return 1
      ;;
  esac
}

check_optional_conda_library() {
  local library_pattern="$1"

  if printf '%s\n' "${ldd_output}" | grep -Eq "[[:space:]]${library_pattern}[^[:space:]]*[[:space:]]+=>"; then
    check_conda_library "${library_pattern}"
  else
    echo "${library_pattern} is not dynamically linked (likely removed by --as-needed)"
  fi
}

check_conda_library 'libstdc\+\+\.so\.6'
check_conda_library 'libgcc_s\.so\.1'
check_conda_library 'libz\.so'
check_conda_library 'libcurl\.so'
check_conda_library 'libhdf5\.so'
check_conda_library 'libgomp\.so\.1'
check_optional_conda_library 'libtbb\.so'
