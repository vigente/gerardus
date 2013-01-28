# Distributed under the Boost Software License, Version 1.0.
# See http://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 2.8.8)
project(BoostMath CXX)

include(Ryppl)
include(RypplPrecompileHeader)

ryppl_find_and_use_package(BoostConfig)

include_directories(
  include
  src/tr1
  )

set(C99_SOURCES
  acosh
  asinh
  atanh
  cbrt
  copysign
  erfc
  erf
  expm1
  fmax
  fmin
  fpclassify
  hypot
  lgamma
  llround
  log1p
  lround
  nextafter
  nexttoward
  round
  tgamma
  trunc
  )

set(TR1_SOURCES
  assoc_laguerre
  assoc_legendre
  beta
  comp_ellint_1
  comp_ellint_2
  comp_ellint_3
  cyl_bessel_i
  cyl_bessel_j
  cyl_bessel_k
  cyl_neumann
  ellint_1
  ellint_2
  ellint_3
  expint
  hermite
  laguerre
  legendre
  riemann_zeta
  sph_bessel
  sph_legendre
  sph_neumann
  )

function(add_math_library name suffix)
  set(sources)
  foreach(src ${ARGN})
    list(APPEND sources "src/tr1/${src}${suffix}.cpp")
  endforeach(src)
  ryppl_add_pch(${name} sources
    <boost/math/special_functions.hpp>
    )
  add_library(boost_math_${name}
    ${sources}
    )
  ryppl_add_pch_to_target(boost_math_${name}
    ${PCH_HEADER}
    )
endfunction(add_math_library)

add_math_library(tr1  ""  ${TR1_SOURCES})
add_math_library(tr1f "f" ${TR1_SOURCES})
add_math_library(tr1l "l" ${TR1_SOURCES})
add_math_library(c99  ""  ${C99_SOURCES})
add_math_library(c99f "f" ${C99_SOURCES})
add_math_library(c99l "l" ${C99_SOURCES})

# ryppl_export(
#   TARGETS
#     boost_math_tr1
#     boost_math_tr1f
#     boost_math_tr1l
#     boost_math_c99
#     boost_math_c99f
#     boost_math_c99l
#   DEPENDS
#     BoostAlgorithm
#     BoostArray
#     BoostConceptCheck
#     BoostConfig
#     BoostConversion
#     BoostCore
#     BoostException
#     BoostFormat
#     BoostFunction
#     BoostFusion
#     BoostInteger
#     BoostMPL
#     BoostRegex
#     BoostSmartPtr
#     BoostStaticAssert
#     BoostTR1
#     BoostTuple
#     BoostTypeTraits
#     BoostUBLAS
#     BoostUtility
#   RECOMMENDS
#     BoostTest
#   INCLUDE_DIRECTORIES
#     include
#   )
