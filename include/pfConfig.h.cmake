/*
 *  Thierry Schuepbach
 *  Vital-IT, Swiss Institute of Bioinformatics
 *
 *  Copyright (C) 2011-2017 Thierry Schuepbach
 */
#ifndef _PFCONFIG_H
#define _PFCONFIG_H

/*
 ********************************************************************************
 * PfTools versioning
 ********************************************************************************
 */
//#include "pfVersion.h"
/*
 ********************************************************************************
 * PfTools integer format
 ********************************************************************************
 */
#cmakedefine USE_32BIT_INTEGER
#ifdef USE_32BIT_INTEGER
# define __USE_32BIT_INTEGER__
#endif

/*
 ********************************************************************************
 * Function choice enumerator
 ********************************************************************************
 */
enum Version { SSE2=0, SSE41=1};

 /*
 ********************************************************************************
 * Enforce inline function
 ********************************************************************************
 */
#ifndef __ALWAYS_INLINE
# ifdef __GNUC__
# ifdef __clang__
#    define  __ALWAYS_INLINE __attribute__((__gnu_inline__, __always_inline__))
#  else
#    define  __ALWAYS_INLINE __attribute__((__gnu_inline__, __always_inline__, __artificial__))
#  endif
#  define __inline inline
# else
#  define  __ALWAYS_INLINE __attribute__((__always_inline__))
# endif
#endif

 /*
 ********************************************************************************
 * Affinity
 ********************************************************************************
 */
#cmakedefine PRF_USE_AFFINITY

/*
 ********************************************************************************
 * Core functionalities
 ********************************************************************************
 */
#cmakedefine PRF_CORE_PCRE
#cmakedefine PRF_CORE_HEURISTIC
#cmakedefine PRF_CORE_STD
#cmakedefine PRF_CORE_REPEAT
#cmakedefine PRF_CORE_FPGA
#cmakedefine PRF_CORE_ZONE
#cmakedefine PRF_CORE_EXT_PROFILE

/*
 ********************************************************************************
 * Input format
 ********************************************************************************
 */
#cmakedefine PRF_INPUT_FASTA
#cmakedefine PRF_INPUT_HDF5
#cmakedefine PRF_INPUT_PBBAM

/*
 ********************************************************************************
 * Specific output format
 ********************************************************************************
 */
#cmakedefine PRF_OUTPUT_PDF
#cmakedefine PRF_OUTPUT_GRAPHICS
#cmakedefine PRF_OUTPUT_DATA

#cmakedefine PRF_OUTPUT_FORMAT_FASTA
#cmakedefine PRF_OUTPUT_FORMAT_INTERPRO
#cmakedefine PRF_OUTPUT_FORMAT_ONELINE
#cmakedefine PRF_OUTPUT_FORMAT_XPSA 
#cmakedefine PRF_OUTPUT_FORMAT_TSV 
#cmakedefine PRF_OUTPUT_FORMAT_FASEARCH 
#cmakedefine PRF_OUTPUT_FORMAT_INCMATCH 
#cmakedefine PRF_OUTPUT_FORMAT_PFSCAN
#cmakedefine PRF_OUTPUT_FORMAT_PSMAKER
#cmakedefine PRF_OUTPUT_FORMAT_SIMPLE
#cmakedefine PRF_OUTPUT_FORMAT_CLASSIFICATION

/*
 ********************************************************************************
 * Mapping database to memory or opening it through libC
 ********************************************************************************
 */
#cmakedefine PRF_USE_MMAP
#ifndef PRF_USE_MMAP
# define SETUP_DATABASE_ACCESS(fileoraddress) FILE* inSequence = fopen((fileoraddress), "r")
# define GET_DATABASE_SEQUENCE(dest, offsets, id) ReadSequenceIndex((dest), (size_t) (id), inSequence, (offsets))
# define UNSET_DATABASE_ACCESS() fclose(inSequence)
#else
# include <sys/mman.h>
# define SETUP_DATABASE_ACCESS(fileoraddress) const char * const restrict SequenceFileMap = (fileoraddress);
# ifndef MMAP_DEBUG
#  define GET_DATABASE_SEQUENCE(dest, offsets, id) MMAP_ReadSequenceIndex((dest), (size_t) (id), SequenceFileMap, (offsets), 0)
# else
#  define GET_DATABASE_SEQUENCE(dest, offsets, id) MMAP_ReadSequenceIndex((dest), (size_t) (id), SequenceFileMap, (offsets), 0\
                                                   ,((struct ThreadData*) _Data)->threadId, 0, *(((struct ThreadData*) _Data)->maplength))
# endif
# define UNSET_DATABASE_ACCESS()
#endif

/*
 ********************************************************************************
 * Import/Export definition
 ********************************************************************************
 */
#ifdef USINGDLL
  #if defined ( WIN32 )
// Visual C/C++, Borland, MinGW and Watcom
    #if defined ( __VISUALC__ ) || defined ( _MSC_VER ) || defined ( __BORLANDC__ ) || defined ( __GNUC__ ) || defined ( __WATCOMC__ )
      #define PFEXPORT    __declspec( dllexport )
      #define PFIMPORT    __declspec( dllimport )
    #else
      #define PFEXPORT
      #define PFIMPORT
    #endif
  #elif defined ( __CYGWIN__ )
    #define PFEXPORT    __declspec( dllexport )
    #define PFIMPORT    __declspec( dllimport )
  #elif defined ( __GNUC__ ) && __GNUC__ > 3
// Follow ideas in http://gcc.gnu.org/wiki/Visibility for GCC version 4.x
// The following forces exported symbols specifically designated with
// PFEXPORT to be visible.
    #define PFEXPORT    __attribute__ ( ( visibility( "default" ) ) )
    #define PFIMPORT
  #endif
#endif

// For an unknown compiler or static built we clear the macros
#ifndef PFEXPORT
  #define PFEXPORT
  #define PFIMPORT
#endif

// The IMPEXP macros will always be set to DLLIMPORT (even for
// the static library, but DLLIMPORT is empty in this case), if
// cmake didn't set the corresponding macro xxxx_EXPORTS when the
// corresponding library is built (DLLIMPEXP is set to DLLEXPORT
// then)
#if defined ( pftools_EXPORTS )
  #define PFIMPEXP    PFEXPORT
  #define PFIMPEXP_DATA( type )    PFEXPORT type
#else
  #define PFIMPEXP    PFIMPORT
  #define PFIMPEXP_DATA( type )    PFIMPORT type
#endif


#endif /* _PFCONFIG_H */
