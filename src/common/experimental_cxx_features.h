#pragma once

#include "common/config.h"

#ifdef HAVE_optional
#include <optional>
#elif defined(HAVE_EXPERIMENTAL_optional)
#include <experimental/optional>
#else
#error "This program needs <optional>"
#endif

#ifdef HAVE_string_view
#include <string_view>
#elif defined(HAVE_EXPERIMENTAL_string_view)
#include <experimental/string_view>
#else
#error "This program needs <string_view>"
#endif

#ifdef HAVE_filesystem
#include <filesystem>
#elif defined(HAVE_EXPERIMENTAL_filesystem)
#include <experimental/filesystem>
#else
#error "This program needs <filesystem>"
#endif

#include <array>
#ifdef HAVE_EXPERIMENTAL_array
#include <experimental/array>
namespace std
{
using experimental::make_array;
}
#endif