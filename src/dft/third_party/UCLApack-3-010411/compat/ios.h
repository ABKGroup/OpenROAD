#pragma once

#include <ios>

// Compatibility shim for legacy iostream headers that used `ios::in`,
// `ios::out`, `ios::beg`, etc. (pre-standard <fstream.h>/<iostream.h> style).
//
// This is intentionally minimal; only the flags used by the UCLApack sources
// are provided.
struct ios
{
  using openmode = std::ios_base::openmode;
  using seekdir = std::ios_base::seekdir;

  static constexpr openmode in = std::ios_base::in;
  static constexpr openmode out = std::ios_base::out;
  static constexpr openmode app = std::ios_base::app;
  static constexpr openmode ate = std::ios_base::ate;
  static constexpr openmode trunc = std::ios_base::trunc;
  static constexpr openmode binary = std::ios_base::binary;

  static constexpr seekdir beg = std::ios_base::beg;
  static constexpr seekdir cur = std::ios_base::cur;
  static constexpr seekdir end = std::ios_base::end;

  // Deprecated legacy flag; for reads, `ios::in` already implies "no create".
  static constexpr openmode nocreate = openmode(0);
};

