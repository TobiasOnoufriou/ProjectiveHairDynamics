// Copyright 2022 Tobias Onoufriou
// Eigen configuration shared by PD-Core and everything that links it.
//
// Force-included (premake `forceincludes`) ahead of Eigen in every source of
// PD-Core, the sandbox, the tests and the tools. It must be identical in all of
// them: Eigen's templates are instantiated in both the library and its users,
// and differing eigen_assert definitions would break the one-definition rule.
//
// eigen_assert only evaluates its condition while pd::g_eigen_asserts is set
// (on by default in Debug builds, off in Release). With
// EIGEN_RUNTIME_NO_MALLOC a test can therefore switch Eigen's assertions, and
// with Eigen::internal::set_is_malloc_allowed(false) its heap-allocation check,
// on for a block of code (StepDoesNotAllocate) without slowing anything else.
#pragma once
#ifdef __cplusplus
#include <cstdio>
#include <cstdlib>

namespace pd {
#ifdef _DEBUG
inline bool g_eigen_asserts = true;
#else
inline bool g_eigen_asserts = false;
#endif
}  // namespace pd

#define EIGEN_RUNTIME_NO_MALLOC
#define eigen_assert(x)                                                     \
  do {                                                                      \
    if (pd::g_eigen_asserts && !(x)) {                                      \
      std::fprintf(stderr, "Eigen assertion failed: %s (%s:%d)\n", #x,      \
                   __FILE__, __LINE__);                                     \
      std::abort();                                                         \
    }                                                                       \
  } while (false)
#endif  // __cplusplus
