####################################
#         Helper cmake vars        #
####################################

# Store a lowercase version of CMAKE_BUILD_TYPE. SCREAM's cmake logic should *always* use this,
# to avoid issues related to case (e.g., Debug vs DEBUG)
set (SCREAM_BUILD_TYPE "" CACHE INTERNAL "A case-insensitive (lowercase) variable accessible everywhere in the project to query the build type.")
string(TOLOWER "${CMAKE_BUILD_TYPE}" SCREAM_BUILD_TYPE)

get_property(IS_EKAT_KOKKOS_BUILT GLOBAL PROPERTY EKAT_KOKKOS_BUILT SET)
if (NOT IS_EKAT_KOKKOS_BUILT)
  message (FATAL_ERROR "Error! You must configure Kokkos *before* including this script.\n")
endif ()

# Determine if this is a Cuda build.
string(FIND "${KOKKOS_GMAKE_DEVICES}" "Cuda" cuda_str_pos)
if (${cuda_str_pos} GREATER -1)
  set(CUDA_BUILD TRUE CACHE INTERNAL "")
else ()
  set(CUDA_BUILD FALSE CACHE INTERNAL "")
endif ()

##################################
#         SCREAM defaults        #
##################################

# Compute reasonable defaults. This needs to happen before the CACHE variables are set.

set (DEFAULT_MIMIC_GPU        FALSE)
set (DEFAULT_POSSIBLY_NO_PACK FALSE)

if (CIME_BUILD)
  set(DEFAULT_LIB_ONLY TRUE)
else ()
  set (DEFAULT_LIB_ONLY FALSE)
endif ()

set (DEFAULT_MAX_THREADS  16)
set (DEFAULT_PACK_SIZE 16)

if (CUDA_BUILD)
  # On the GPU, the pack size must be 1
  set(DEFAULT_PACK_SIZE 1)
  set(DEFAULT_MAX_THREADS 1)
else()
  if (SCREAM_BUILD_TYPE STREQUAL "debug")
    set(DEFAULT_MIMIC_GPU TRUE)
  endif()
endif ()

set (DEFAULT_FPMODEL "precise")
if (SCREAM_DEFAULT_MPIRUN_EXE)
  set(DEFAULT_MPIRUN_EXE ${SCREAM_DEFAULT_MPIRUN_EXE})
else()
  set(DEFAULT_MPIRUN_EXE "mpiexec")
endif()

# For some routines, SKX may have better performance with pksize=1
if ("${KOKKOS_GMAKE_ARCH}" STREQUAL "SKX")
  set(DEFAULT_POSSIBLY_NO_PACK TRUE)
endif ()

if (SCREAM_BUILD_TYPE STREQUAL "debug")
  set(DEFAULT_FPMODEL "strict")
endif()

##################################
#         SCREAM options         #
##################################

# Whether SCREAM should use double precision real values
set(SCREAM_DOUBLE_PRECISION TRUE CACHE BOOL "Set to double precision (default True)")

# Whether SCREAM should try to mimic on CPU the way parallel work is distributed on GPU
set(SCREAM_MIMIC_GPU ${DEFAULT_MIMIC_GPU} CACHE BOOL "Mimic GPU to correctness-test inter-column parallelism on non-GPU platform")

# Where scream puts files that are generated for testing purposes
set(SCREAM_TEST_DATA_DIR ${CMAKE_CURRENT_BINARY_DIR}/data CACHE FILEPATH "Location of data files generated by tests")

# Whether scream should error out when an MPI error occurs
set(SCREAM_MPI_ERRORS_ARE_FATAL TRUE CACHE BOOL "Whether MPI errors should abort (default TRUE). If false, errors should be handled.")

# The name of the executable used to launche mpi executables
set(SCREAM_MPIRUN_EXE ${DEFAULT_MPIRUN_EXE} CACHE STRING "The executable name for mpirun")

# Whether SCREAM should build executables or not
set(SCREAM_LIB_ONLY ${DEFAULT_LIB_ONLY} CACHE BOOL "Only build libraries, no exes")

# Whether the pack object should check bounds when accessing entries
set(SCREAM_PACK_CHECK_BOUNDS FALSE CACHE BOOL "If defined, scream::pack objects check indices against bounds")

# The master pack size
set(SCREAM_PACK_SIZE ${DEFAULT_PACK_SIZE} CACHE STRING
  "The number of scalars in a scream::pack::Pack and Mask. Larger packs have good performance on conditional-free loops due to improved caching.")

# With the master pack size determined, we have constraints on the others.
set(SCREAM_SMALL_PACK_SIZE ${SCREAM_PACK_SIZE} CACHE STRING
  "The number of scalars in a scream::pack::SmallPack and SmallMask. Smaller packs can have better performance in loops with conditionals since more of the packs will have masks with uniform value.")
set (SCREAM_POSSIBLY_NO_PACK ${DEFAULT_POSSIBLY_NO_PACK} CACHE BOOL
  "Set possibly-no-pack to this value. You can set it to something else to restore packs on SKX for testing.")

# Scream test configuration options
set(SCREAM_TEST_MAX_THREADS ${DEFAULT_MAX_THREADS} CACHE STRING "Upper limit on threads for threaded tests")
set(SCREAM_TEST_THREAD_INC 1 CACHE STRING "Thread count increment for threaded tests")

# The name of the library providing the dycore for dynamics
set(SCREAM_DYNAMICS_DYCORE "NONE" CACHE STRING "The name of the dycore to be used for dynamics. If NONE, then any code/test requiring dynamics is disabled.")

set(SCREAM_FPMODEL ${DEFAULT_FPMODEL} CACHE STRING "Compiler floating point model")

# The following are options that are either avaiable and processed by CMake *only if* certain conditions
# on other options are met, OR whose default value depends on other options.
# The syntax of cmake_dependent_option is
#  cmake_dependent_option (opt_name "description" default_if_true "conditions" value_if_false)
# where conditions is a semi-colon separated list of boolean expressions, like "VAR1; NOT VAR2; VAR3"
# Unfortunately, there is no cmake_dependent_cache_var, so non-boolean dependent config options need 
# to be manually processed with if/else statements

include (CMakeDependentOption)

# Assuming SCREAM_LIB_ONLY is FALSE (else, no exec is built at all), we provide the option
# of building only baseline-related execs. By default, this option is off (menaing "build everything).
# However, when generating baselines, this can be useful to reduce the amount of stuff compiled.
cmake_dependent_option (SCREAM_BASELINES_ONLY
                        "Whether building only baselines-related executables" OFF
                        "NOT SCREAM_LIB_ONLY" OFF)

if (NOT CUDA_BUILD AND ${SCREAM_PACK_SIZE} EQUAL 1 AND SCREAM_BUILD_TYPE STREQUAL "debug")
  set (DEFAULT_FPE ON)
endif ()
option (SCREAM_FPE "Enable floating point error exception" ${DEFAULT_FPE})
