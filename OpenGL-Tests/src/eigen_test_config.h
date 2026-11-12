// Copyright 2022 Tobias Onoufriou
// Force-included before every source of the test project (premake
// forceincludes), ahead of Eigen.
//
// Release builds compile eigen_assert away, which would also silence the
// EIGEN_RUNTIME_NO_MALLOC check. This eigen_assert only evaluates its condition
// while g_eigen_asserts is set, so a test can switch Eigen's assertions (and
// with Eigen::internal::set_is_malloc_allowed(false) its heap-allocation check)
// on for a block of code without slowing down the other tests and benchmarks.
#pragma once
#include <cstdio>
#include <cstdlib>

namespace pd_test {
inline bool g_eigen_asserts = false;
}  // namespace pd_test

#define EIGEN_RUNTIME_NO_MALLOC
#define eigen_assert(x)                                                     \
  do {                                                                      \
    if (pd_test::g_eigen_asserts && !(x)) {                                 \
      std::fprintf(stderr, "Eigen assertion failed: %s (%s:%d)\n", #x,      \
                   __FILE__, __LINE__);                                     \
      std::abort();                                                         \
    }                                                                       \
  } while (false)
