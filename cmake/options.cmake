OPTION(BUILD_SHARED_LIBS "Enable shared library" ON)
OPTION(BUILD_STATIC_LIBS "Enable static library" OFF)
OPTION(BUILD_TOOLS "Enable tools" OFF)
OPTION(TOOLS_STATIC_BUILD  "This will compile tools in static mode." OFF)

OPTION(PRF_USE_AFFINITY "CPU core affinity will be used for threads." ON)

OPTION(PRF_USE_MMAP "DB files wil be mapped to memory for speed." ON)
OPTION(PRF_USE_32BIT_INTEGER "Builds using 32bit integer format rather than 16 bits." OFF)

IF (BUILD_STATIC_LIBS)
	OPTION(BUILD_STATIC_DEPENDENCIES "Enable entire static library" OFF)
ENDIF(BUILD_STATIC_LIBS)

OPTION(PRF_INPUT_FASTA "Enable FASTA input" ON)
OPTION(PRF_INPUT_HDF5 "Enable PacBio HDF5 input" OFF)
OPTION(PRF_INPUT_PBBAM "Enable PacBio BAM input" OFF)

OPTION(PRF_OUTPUT_PDF "Enable output in pdf file format." OFF)
OPTION(PRF_OUTPUT_GRAPHICS "Enable graphics generation." OFF)
OPTION(PRF_OUTPUT_DATA "Enable direct dump of memory array to file." ON)

OPTION(PRF_OUTPUT_FORMAT_FASTA "Output in fasta format" ON)
OPTION(PRF_OUTPUT_FORMAT_ONELINE "Output in OneLine format" ON)
OPTION(PRF_OUTPUT_FORMAT_XPSA "Output in xPSA format" ON)
OPTION(PRF_OUTPUT_FORMAT_TSV "Output in TSV format" ON)
OPTION(PRF_OUTPUT_FORMAT_INTERPRO "Output in interpro format" OFF)
OPTION(PRF_OUTPUT_FORMAT_FASEARCH "Output in FaSearch format" OFF)
OPTION(PRF_OUTPUT_FORMAT_INCMATCH "Output in IncMatch format" OFF)
OPTION(PRF_OUTPUT_FORMAT_PFSCAN "Output in pfScan format" OFF)
OPTION(PRF_OUTPUT_FORMAT_PSMAKER "Output in psMaker format" OFF)
OPTION(PRF_OUTPUT_FORMAT_SIMPLE "Output in simple format" OFF)
OPTION(PRF_OUTPUT_FORMAT_CLASSIFICATION "Output in Classification format" OFF)
OPTION(PRF_OUTPUT_FORMAT_TEST "Output for test format" OFF)

OPTION(PRF_CORE_MAP "Mapping routines" ON)
OPTION(PRF_CORE_HEURISTIC "Heuristic routines" OFF)
OPTION(PRF_CORE_STD "Standard alignment routines" ON)
OPTION(PRF_CORE_REPEAT "Circular alignment routine" OFF)
OPTION(PRF_CORE_FPGA "Specific FPGA like routine" OFF)
OPTION(PRF_CORE_ZONE "Specific alignment routine to get zones" OFF)
OPTION(PRF_CORE_EXT_PROFILE "Extended profile alignment for genome caller" OFF)

OPTION(PRF_CORE_PCRE "Enable Perl Regex." ON)
OPTION(PRF_CORE_DISPATCH "Enable dispatching to threads." ON)

OPTION(PRF_WRAPPER_JAVA "Addons far Java Native Interface." OFF)
