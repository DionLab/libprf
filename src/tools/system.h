/*******************************************************
                        PFTOOLS
 *******************************************************
  Sep 26, 2011 system.h
 *******************************************************
 (C) 2011 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/
#ifndef SYSTEM_H_
#define SYSTEM_H_
#ifdef __USE_AFFINITY__
#include <sched.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#if defined(_WIN32)
#include <windows.h>
#include <winbase.h>
#else
#include <sys/utsname.h>
#include <sys/types.h>
#include <unistd.h>
#include <pwd.h>
#endif
#ifdef __NUMA__
#  include <numa.h>
#endif


#define MM_MMX    	0x00000001 /* standard MMX */
#define MM_MMXEXT 	0x00000002 /* SSE integer functions or AMD MMX ext */
#define MM_3DNOW  	0x00000004 /* AMD 3DNOW */
#define MM_SSE    	0x00000008 /* SSE functions */
#define MM_SSE2   	0x00000010 /* PIV SSE2 functions */
#define MM_3DNOWEXT  	0x00000020 /* AMD 3DNowExt */
#define MM_SSE3   	0x00000040 /* Prescott SSE3 functions */
#define MM_SSSE3  	0x00000080 /* Conroe SSSE3 functions */
#define MM_SSE41  	0x00000100 /* Penryn SSE41 functions */
#define MM_SSE42  	0x00000200 /* Nehalem SSE42 functions */
#define MM_POPCOUNT 	0x00000400
#define MM_SSE4A	0x00000800 /* AMD Barcelona functions */
#define MM_AVX          0x00001000 /* AVX technology 256 bits */

#ifdef _SYSTEM_TEST
#define __SYSTEM_DEBUG__
#endif

/* TOPOLOGY structure is as follow
 *
 * 2  bits for hyperthreading threads b0, b1 : b1 is hyperthreading only, b0 is non hyperthreading
 * 16 bits for the cores
 * 14 bits for the sockets
 * ------------------------
 * 32 bits
 */


typedef struct SystemInfo {
    /* Operating system data */
    char Username[64];
    char Release[64];
    char Architecture[64];
    char Nodename[64];
    /* Architecture data */
    char CPU_Vendor[13];
    char CPU_Name[48];
    unsigned int Family;
    unsigned int Model;
    unsigned int Extensions;
    /* Topology */
    unsigned int nSockets;
    unsigned int nCores;
    unsigned int nOverallCores;
    _Bool HyperthreadingAvailable;
    _Bool HyperthreadingOn;
    unsigned int nComputeUnit;
    unsigned int nCorePerComputeUnit;

    unsigned int SocketSelectMask;
    unsigned int CoreSelectMask;
    unsigned int HyperthreadingMask;
    unsigned int SocketSelectMaskShift;
    unsigned int HyperthreadingMaskWidth;

    unsigned int * TopologyMasks;
    unsigned int * APICID;
#ifdef __NUMA__
    unsigned int nNodes;
    unsigned int nCpusPerNode;
    _Bool NumaAble;
#endif
} __attribute__((aligned(16))) SystemInfo;


extern __inline unsigned long __attribute__((__gnu_inline__, __always_inline__)) createMask(unsigned int numEntries, unsigned int *maskWidth)
{
    unsigned int i;
    unsigned long k = ((unsigned long) numEntries) * 2L - 1L;
    unsigned long index;
    const unsigned int one = 1;
    unsigned int AllZero;
#if (defined(__x86_64__) || defined(__amd64__))&& !defined(_WIN32)
    __asm__ __volatile__ (
      "xorl   %0, %0; "
      "bsrq  %2,%1; "
      "cmovzl %3,%0;"
      : "=r"(AllZero), "=r"(index)
      : "r"(k), "m"(one)
    );
#else
 __asm__ __volatile__ (
      "xorl   %0, %0; "
      "bsrl  %2,%1; "
      "cmovzl %3,%0;"
      : "=r"(AllZero), "=r"(index)
      : "r"(k), "m"(one)
    );
//#error Code still needs to be ported to 32 bit architecture
#endif
    if (AllZero) {
      *maskWidth = 0;
      return 0L;
    }
    *maskWidth = (unsigned int) index;
    if (index == 63L) {
      return -1L;
    } else {
      return (1L << index) - 1L;
    }
}


extern __inline void __attribute__((__gnu_inline__, __always_inline__)) freeSystemInfo(SystemInfo * const info)
{
  if (info->APICID) {
    free(info->APICID);
  }
}

#if defined(__USE_AFFINITY__) && !defined(_WIN32)
//extern __inline unsigned int __attribute__((__gnu_inline__, __always_inline__))
unsigned int
getMasks(SystemInfo * const info, const unsigned int SocketId, const unsigned int CoreId, const unsigned int ThreadId, cpu_set_t * * const Masks)
{
  unsigned int count = 0;
  unsigned int ValidMask;

  if (SocketId == 0 || CoreId == 0 || ThreadId == 0 || info->TopologyMasks == NULL) return 0;
  ValidMask  = ThreadId & 0x3;
  ValidMask |= (CoreId & 0xFFFF) << 2;
  ValidMask |= (SocketId & 0x3FFF) << 18;
#ifdef __SYSTEM_DEBUG__
  fputs("Valid Mask  : ", stdout);
  for (unsigned int j=0; j<32; ++j) {
    if (ValidMask & (1 << (31 - j)) )
      fputs("1",stdout);
    else
      fputs("0", stdout);
    if (j == 13 || j == 29) fputs(" ",stdout);
  }
  fputs("\n",stdout);
#endif
  for (unsigned int thread=0; thread<info->nOverallCores; ++thread) {
    register const unsigned int tmp = info->TopologyMasks[thread] & ValidMask;
    if ((tmp & 0x3) && (tmp & 0x3FFFC) && (tmp & 0xFFFC0000)) ++count;
  }

  if (count == 0) return 0;

  *Masks = (cpu_set_t*) malloc(count*sizeof(cpu_set_t));
  if (*Masks == NULL) return 0;
  count = 0;
  for (unsigned int thread=0; thread<info->nOverallCores; ++thread) {
    register const unsigned int tmp = info->TopologyMasks[thread] & ValidMask;
    if ((tmp & 0x3) && (tmp & 0x3FFFC) && (tmp & 0xFFFC0000)) {
       CPU_ZERO(*Masks + count);
       CPU_SET(thread, *Masks + count);
#ifdef __SYSTEM_DEBUG__
       printf("Thread %2u   : ", thread);
       for (unsigned int j=0; j<32; ++j) {
	if (info->TopologyMasks[thread] & (1 << (31 - j)) )
	  fputs("1",stdout);
	else
	  fputs("0", stdout);
	if (j == 13 || j == 29) fputs(" ",stdout);
      }
      fputs("\n",stdout);
#endif
      ++count;
    }
  }
  return count;
}
#endif
//extern __inline void __attribute__((__gnu_inline__, __always_inline__)) 
void
getSystemInfo(SystemInfo * const info)
{
#if defined(_WIN32)

#else
  struct utsname uts_name;
#endif
  char username[16] __attribute__((aligned(16)));
  time_t creation_time;

  /* clear all data */
  memset(info, 0, sizeof(SystemInfo));

  /*
   *  OPERATING SYSTEM INFORMATIONS ------------------------------------------------
   */
#if defined(_WIN32)

#else
  uname(&uts_name);
  size_t tmp_size = strlen(uts_name.release);
	if (tmp_size<63) {
		strncpy(info->Release, uts_name.release, 64);
  } else {
		memcpy(info->Release, uts_name.release, 60);
		info->Release[60] = '.';
		info->Release[61] = '.';
		info->Release[62] = '.';
		info->Nodename[63] = '\0';
  }
#ifndef NO_USERNAME
  uid_t uid = getuid();
  struct passwd *pass= getpwuid(uid);
  if ( pass != NULL) {
      if (strlen(pass->pw_gecos) > 0) {
				strncpy(info->Username,pass->pw_gecos,64);
      } else {
				strncpy(info->Username,pass->pw_name,64);
      }
  }
#else
  const char text[] = "Static linking prevent username query";
  strcpy(info->Username, text);
#endif

  tmp_size = strlen(uts_name.nodename);
  if (tmp_size<63) {
		strncpy(info->Nodename, uts_name.nodename, 64);
  } else {
		memcpy(info->Nodename, uts_name.nodename, 60);
		info->Nodename[60] = '.';
		info->Nodename[61] = '.';
		info->Nodename[62] = '.';
		info->Nodename[63] = '\0';
  }

  tmp_size = strlen(uts_name.machine);
  if (tmp_size<63) {
      strncpy(info->Architecture, uts_name.machine, 64);
  } else {
		memcpy(info->Architecture, uts_name.machine, 60);
		info->Architecture[60] = '.';
		info->Architecture[61] = '.';
		info->Architecture[62] = '.';
		info->Nodename[63] = '\0';
  }
  for (size_t i=0; i<13; ++i) info->CPU_Vendor[i] = '\0';

  /* Available number of logical processor seen by OS */
  info->nOverallCores = (unsigned int) sysconf(_SC_NPROCESSORS_CONF);
#endif
  /*
   * ARCHITECTURE INFORMATION -----------------------------------------------
   */

  const char NOBrand[] = "Not supported by this cpu";
  strcpy(info->CPU_Name, NOBrand);

  unsigned long a,c;  
  __asm__ __volatile__ (
    /* See if CPUID instruction is supported ... */
    /* ... Get copies of EFLAGS into eax and ecx */
    "pushf\n\t"
    "popq %%rax\n\t"
    "movq %0, %1\n\t"

    /* ... Toggle the ID bit in one copy and store */
    /*     to the EFLAGS reg */
    "xorq $0x200000, %0\n\t"
    "pushq %%rax\n\t"
    "popf\n\t"

    /* ... Get the (hopefully modified) EFLAGS */
    "pushf\n\t"
    "popq %%rax\n\t"
    : "=a" (a), "=c" (c)
    :
    : "cc"
    );

  unsigned int rval = 0;
  if (a != c) {
    int eax_value, ebx_value, ecx_value, edx_value;

    /* Get the CPU vendor string */
		__asm__ __volatile__ (
			"xorl %%eax, %%eax;"
	    "xorl %%ecx, %%ecx;"
	    "cpuid ;"
	    "movl  %%ebx,   (%2);"
	    "movl  %%edx,  4(%2);"
	    "movl  %%ecx,  8(%2);"
	    //"mov  %%ebx, 4(%0) ;" 
	    : "=&a"(eax_value)
	    : "0" (eax_value), "r"(info->CPU_Vendor)
	    : "memory", /*"%ebp",*/ "%ebx", "%ecx", "%edx");
		
    /* if available query SSE and MMX functionalities */
    if (eax_value >= 1) {
      eax_value = 1;
      __asm__ __volatile__ (
	      "cpuid \n\t"
	      : "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
	      : "0"(eax_value)
	      );

      if (edx_value & (1<<23))    rval |= MM_MMX;
      if (edx_value & (1<<25))    rval |= MM_MMXEXT | MM_SSE;
      if (edx_value & (1<<26))    rval |= MM_SSE2;
      if (ecx_value & 1)	  rval |= MM_SSE3;
      if (ecx_value & (1<<9) )    rval |= MM_SSSE3;
      if (ecx_value & (1<<19))    rval |= MM_SSE41;
      if (ecx_value & (1<<20))    rval |= MM_SSE42;
      if (ecx_value & (1<<23))    rval |= MM_POPCOUNT;
      if (ecx_value & (1<<28))    rval |= MM_AVX;

      info->Family = (( eax_value >> 8 ) & 0xF) + (( eax_value >> 20 ) & 0xFF);
      info->Model = (( eax_value >> 4 ) & 0xF) || ((( eax_value >> 16 ) & 0xF) << 4);
    }

     __asm__ __volatile__ (
	      "xorl %%ecx, %%ecx\n\t "
	      "cpuid \n\t"
	      : "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
	      : "0"(0x80000000)
	      );
    const unsigned int MaxLeaf = eax_value;
#ifdef __SYSTEM_DEBUG__
    fprintf(stderr, "Maximum extended leaf : 0x%8xh\n", MaxLeaf);
#endif
    if (MaxLeaf >= 0x80000001) {
      __asm__ __volatile__ (
	      "cpuid \n\t"
	      : "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
	      : "0"(0x80000001)
	      );

      if (edx_value & (1<<31))  rval |= MM_3DNOW;
      if (edx_value & (1<<30))  rval |= MM_3DNOWEXT;
      if (edx_value & (1<<23))  rval |= MM_MMX;
      if (edx_value & (1<<22))  rval |= MM_MMXEXT;
      if (ecx_value & (1<<6))   rval |= MM_SSE4A;
    }

    if (MaxLeaf >= 0x80000004) {
      /* Get the CPU brand string */
      __asm__ __volatile__ (
	    "xorl  %%ecx, %%ecx\n\t"
	    "movl  $0x80000002, %%eax\n\t"
	    "cpuid \n\t"
	    "movl  %%eax,   (%0)\n\t"
	    "movl  %%ebx,  4(%0)\n\t"
	    "movl  %%ecx,  8(%0)\n\t"
	    "movl  %%edx, 12(%0)\n\t"
	    "movl  $0x80000003, %%eax\n\t"
	    "cpuid \n\t"
	    "movl  %%eax, 16(%0)\n\t"
	    "movl  %%ebx, 20(%0)\n\t"
	    "movl  %%ecx, 24(%0)\n\t"
	    "movl  %%edx, 28(%0)\n\t"
	    "movl  $0x80000004, %%eax\n\t"
	    "cpuid \n\t"
	    "movl  %%eax, 32(%0)\n\t"
	    "movl  %%ebx, 36(%0)\n\t"
	    "movl  %%ecx, 40(%0)\n\t"
	    "movl  %%edx, 44(%0)\n\t"
	    :
	    : "r"(info->CPU_Name)
	    : "memory", "%eax", "%ebx", "%ecx", "%edx");
      {
				/* Correct space at the beginning */
				char * ptr = info->CPU_Name;
				size_t i = 0;
				while (ptr[i] == ' ' && i < 48 ) ++i;
				while (i < 48 ) { *ptr++ = info->CPU_Name[i++]; }
      }
    }

    info->Extensions = rval;

    /*
     * TOPOLOGY INFORMATIONS ---------------------------------------------------------
     */

    /* Intel Topology*/
    if (info->CPU_Vendor[0] == 'G' && info->CPU_Vendor[1] == 'e')
    {
      /* Do we support Leaf 11 (0xB) */
      _Bool UseLeafB = false;
      __asm__ __volatile__ (
	      "xorl %%eax, %%eax\n\t"
	      "xorl %%ecx, %%ecx\n\t"
	      "cpuid \n\t"
	      : "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
	      );
      const unsigned int MaxCPUID = eax_value;
      if (eax_value >= 0xB) {
				eax_value = 0xB;
				__asm__ __volatile__ (
	      "xorl %%ecx, %%ecx\n\t"
	      "cpuid \n\t"
	      : "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
	      : "0"(eax_value), "2"(ecx_value)
	      );
				UseLeafB = (ebx_value != 0);
      }
      /* Use HWMT hyperthreading feature flag to treat different configurations */
      eax_value = 1;
      __asm__ __volatile__ (
	      "xorl %%ecx, %%ecx\n\t"
	      "cpuid \n\t"
	      : "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
	      : "0"(eax_value)
	      );

      if (edx_value & (1<<28)) {
				int wasCoreReported = 0;
				int wasThreadReported = 0;
				unsigned int ThreadPerCore = 1;
				unsigned int MaxCore = 1;

				if (UseLeafB) {
					int subLeaf = 0, levelType, levelShift;
					unsigned int coreplusSMT_Mask;
					do {
						eax_value = 0xB;
						ecx_value = subLeaf;
						__asm__ __volatile__ (
							"cpuid \n\t"
							: "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
							: "0"(eax_value), "2"(ecx_value)
							);
						if (ebx_value == 0) break;

						levelType  = (ecx_value >> 8) & 0xFF;
						levelShift = (eax_value & 0xF);
						switch (levelType)
						{
							case 1:
					//level type is SMT, so levelShift is the SMT_Mask_Width
					info->HyperthreadingMask = ~((-1) << levelShift);
					info->HyperthreadingMaskWidth = levelShift;
					wasThreadReported = 1;
					break;
							case 2: //level type is Core, so levelShift is the CorePlsuSMT_Mask_Width
					coreplusSMT_Mask = ~((-1) << levelShift);
					info->SocketSelectMaskShift = levelShift;
					info->SocketSelectMask = (-1) ^ coreplusSMT_Mask;
					wasCoreReported = 1;
					break;
							default:
					// handle in the future
					break;
						}
						++subLeaf;
					} while (1);

					if (wasThreadReported && wasCoreReported)
					{
						info->CoreSelectMask = coreplusSMT_Mask ^ info->HyperthreadingMask;
					}
					else if (!wasCoreReported && wasThreadReported)
					{
						info->CoreSelectMask = 0;
						info->SocketSelectMaskShift = info->HyperthreadingMaskWidth;
						info->SocketSelectMask = (-1) ^ info->HyperthreadingMask;
					}
					else //(case where !wasThreadReported)
					{
						// throw an error, this should not happen if hardware function normally
						fputs("Error in hardware info extraction\n", stderr);
						exit(1);
					}
					info->HyperthreadingAvailable = true;
					info->HyperthreadingOn = (wasThreadReported > 0);

				} else {
					const unsigned int MaxCorePlusThread = (ebx_value >> 16) & 0xFF;
					if (MaxCPUID >= 4) {
						eax_value = 4;
						__asm__ __volatile__ (
							"xorl %%ecx, %%ecx\n\t"
							"cpuid \n\t"
							: "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
							: "0"(eax_value)
							);
						MaxCore = 1 + ( ( eax_value >> 26 ) & 0x3F );
						ThreadPerCore = MaxCorePlusThread / MaxCore;
						info->HyperthreadingAvailable = true;
						info->HyperthreadingOn = (MaxCorePlusThread > MaxCore ) ? true : false;

					} else {
						MaxCore = 1;
						ThreadPerCore = MaxCorePlusThread;

						/* Check whether BIOS is preventing Hyperthreading */
						eax_value = 0x80000000;
						__asm__ __volatile__ (
							"xorl %%ecx, %%ecx\n\t"
							"cpuid \n\t"
							: "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
							: "0"(eax_value)
							);
						if ( MaxCPUID <= 4 && eax_value > 0x80000004) {
							info->HyperthreadingAvailable = true;
						}
					}

					/* Create the masks */
					info->HyperthreadingMask = createMask(ThreadPerCore, &(info->HyperthreadingMaskWidth));
					info->CoreSelectMask = createMask(MaxCore, &(info->SocketSelectMaskShift));
					info->SocketSelectMaskShift += info->HyperthreadingMaskWidth;
					info->CoreSelectMask <<= info->HyperthreadingMaskWidth;
					info->SocketSelectMask = (-1) ^ (info->CoreSelectMask | info->HyperthreadingMask);
				}
      } else {
				/* Prior to Hyperthreading Technology only one thread per core */
				info->SocketSelectMask = -1;
      }
#if defined(__USE_AFFINITY__) && !defined(_WIN32)
      /* Retrieve APIC data for each thread */
      cpu_set_t Mask;
      cpu_set_t BackupMask;
      if ( sched_getaffinity(0, sizeof(cpu_set_t), &BackupMask)) {
				fputs("Error getting affinity!\n",stderr);
				exit(1);
      }

      info->APICID = (unsigned int*) malloc(info->nOverallCores*(sizeof(unsigned int) + sizeof(unsigned int)));
      info->TopologyMasks = &(info->APICID[info->nOverallCores]);

      if (UseLeafB) {
				for (unsigned int thread =0; thread<info->nOverallCores; ++thread) {
					CPU_ZERO(&Mask);
					CPU_SET(thread, &Mask);
					if ( sched_setaffinity(0, sizeof(cpu_set_t), &Mask)) {
							fputs("Error setting affinity!\n",stderr);
							exit(1);
					}
					__asm__ __volatile__ (
							"xorl %%ecx, %%ecx\n\t"
							"cpuid\n\t"
							: "=d"(edx_value)
							: "a"(0xB)
							: "memory", "ebx", "ecx"
						);
					info->APICID[thread] = edx_value;
				}
      } else {
				for (unsigned int thread =0; thread<info->nOverallCores; ++thread) {
					CPU_ZERO(&Mask);
					CPU_SET(thread, &Mask);
					if ( sched_setaffinity(0, sizeof(cpu_set_t), &Mask)) {
							fputs("Error setting affinity!\n",stderr);
							exit(1);
					}
					__asm__ __volatile__ (
							"xorl %%ecx, %%ecx\n\t"
							"cpuid\n\t"
							: "=b"(ebx_value)
							: "a"(0x1)
							: "memory", "%ecx", "%edx"
						);
					unsigned int tmp = ( ebx_value >> 24 ) & 0xFF;
					info->APICID[thread] = tmp;
				}
      }

      if ( sched_setaffinity(0, sizeof(cpu_set_t), &BackupMask)) {
				fputs("Error setting affinity!\n",stderr);
				exit(1);
      }
#ifdef __SYSTEM_DEBUG__
      printf("HTT Mask (%2u)    : ",info->HyperthreadingMaskWidth);
      for (unsigned int j=0; j<32; ++j) {
				if ( info->HyperthreadingMask & (1 << (31 - j)) )
					fputs("1",stdout);
				else
					fputs("0", stdout);
      }
      fputs("\nCore Mask        : ",stdout);
      for (unsigned int j=0; j<32; ++j) {
				if ( info->CoreSelectMask & (1 << (31 - j)) )
					fputs("1",stdout);
				else
					fputs("0", stdout);
      }
      printf("\nSocket Mask (%2u) : ", info->SocketSelectMaskShift);
      for (unsigned int j=0; j<32; ++j) {
				if ( info->SocketSelectMask & (1 << (31 - j)) )
					fputs("1",stdout);
				else
					fputs("0", stdout);
      }
      fputs("\n",stdout);
#endif
      register unsigned int nSockets = 0;
      register unsigned int nCores = 0;
      for (unsigned int thread =0; thread<info->nOverallCores; ++thread) {
				register unsigned int TMask = 1 << (info->APICID[thread] & info->HyperthreadingMask);
				register unsigned int tmp = (info->APICID[thread] & info->CoreSelectMask) >> info->HyperthreadingMaskWidth;
				if (tmp > nCores) nCores = tmp;
				TMask |= 1 << (2 + tmp);
				tmp = (info->APICID[thread] & info->SocketSelectMask) >> info->SocketSelectMaskShift;
				if (tmp > nSockets) nSockets = tmp;
				TMask |= 1 << ( 18 + tmp);
				info->TopologyMasks[thread] = TMask;
      }
      info->nCores = 1 + nCores;
      info->nSockets = 1 + nSockets;
#endif
    }
    /* AMD Topology*/
    else if (info->CPU_Vendor[0] == 'A' && info->CPU_Vendor[1] == 'u')
    {
      /* Use HTT hyperthreading feature flag to test single or multicore architecture */
      eax_value = 1;
      __asm__ __volatile__ (
	      "xorl %%ecx, %%ecx\n\t"
	      "cpuid \n\t"
	      : "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
	      : "0"(eax_value)
	      );

      if (edx_value & (1<<28)) {
				/* Get the number of logical core per processor */
				info->nCores = 1 + ((ebx_value >> 16) & 0xF);
      }

      /*Get the largest extension leaf supported */
      eax_value = 0x80000000;
      __asm__ __volatile__ (
	    "xorl %%ecx, %%ecx\n\t"
	    "cpuid \n\t"
	    : "=a"(eax_value)
	    : "0"(eax_value)
	    : "%ebx", "%ecx", "%edx"
      );

      /* Bulldozer */
      const _Bool IsBulldozer = eax_value >= 0x8000001E ? true : false;

      /* Get the number of physical core per processor */
      eax_value = 0x80000008;
      __asm__ __volatile__ (
	    "xorl %%ecx, %%ecx\n\t"
	    "cpuid \n\t"
	    : "=a"(eax_value), "=b"(ebx_value), "=c"(ecx_value), "=d"(edx_value)
	    : "0"(eax_value)
      );

      info->nCores = 1 + (ecx_value & 0xF );
      unsigned int coreplusSMT_MaskWidth = (ecx_value >> 12) & 0xF;
#if defined(__USE_AFFINITY__) && !defined(_WIN32)
      /* Retrieve APIC data for each thread */
      cpu_set_t Mask;
      cpu_set_t BackupMask;
      if ( sched_getaffinity(0, sizeof(cpu_set_t), &BackupMask)) {
				fputs("Error getting affinity!\n",stderr);
				exit(1);
      }

      info->APICID = (unsigned int*) malloc(info->nOverallCores*(sizeof(unsigned int) + sizeof(unsigned int)));
      info->TopologyMasks = &(info->APICID[info->nOverallCores]);

      if (!IsBulldozer) {
				for (unsigned int thread =0; thread<info->nOverallCores; ++thread) {
					CPU_ZERO(&Mask);
					CPU_SET(thread, &Mask);
					if ( sched_setaffinity(0, sizeof(cpu_set_t), &Mask)) {
							fputs("Error setting affinity!\n",stderr);
							exit(1);
					}
					__asm__ __volatile__ (
							"xorl %%ecx, %%ecx\n\t"
							"cpuid\n\t"
							: "=b"(ebx_value)
							: "a"(0x1)
							: "%ecx", "%edx"
						);
					unsigned int tmp = ( ebx_value >> 24 ) & 0xFF;
					info->APICID[thread] = tmp;
				}
      } else {
				/* Get number of core per compute unit assuming all are identical */
				__asm__ __volatile__ (
					"xorl %%ecx, %%ecx\n\t"
					"cpuid\n\t"
					: "=b"(ebx_value)
					: "a"(0x8000001E)
					: "%ecx", "%edx"
				);
				info->nCorePerComputeUnit = 1 + ( ( ebx_value >> 8) & 0xFF );
				unsigned int nCU = 0;
				for (unsigned int thread =0; thread<info->nOverallCores; ++thread) {
					CPU_ZERO(&Mask);
					CPU_SET(thread, &Mask);
					if ( sched_setaffinity(0, sizeof(cpu_set_t), &Mask)) {
							fputs("Error setting affinity!\n",stderr);
							exit(1);
					}

					__asm__ __volatile__ (
							"xorl %%ecx, %%ecx\n\t"
							"cpuid\n\t"
							: "=b"(ebx_value)
							: "a"(0x1)
							: "%ecx", "%edx"
						);
					unsigned int tmp = ( ebx_value >> 24 ) & 0xFF;
					info->APICID[thread] = tmp;
#ifdef __SYSTEM_DEBUG__
					fprintf(stderr,"APIC DATA (%2u)   : ", thread);
					for (unsigned int j=0; j<8; ++j) {
						if ( tmp & (1 << (7 - j)) )
							fputs("1",stderr);
						else
							fputs("0", stderr);
					}
#endif
					__asm__ __volatile__ (
							"xorl %%ecx, %%ecx\n\t"
							"cpuid\n\t"
							: "=b"(ebx_value)
							: "a"(0x8000001E)
							: "%ecx", "%edx"
						);

					tmp = ebx_value & 0xFF;
					if (tmp > nCU) nCU = tmp;
#ifdef __SYSTEM_DEBUG__
					fputs(" CU ",stderr);
					for (unsigned int j=0; j<8; ++j) {
						if ( tmp & (1 << (7 - j)) )
							fputs("1",stderr);
						else
							fputs("0", stderr);
					}
					fputs("\n", stderr);
#endif
				}
				info->nComputeUnit = 1 + nCU;
			}

      /* Create the masks */
      if ( info->nCorePerComputeUnit ) {
				info->HyperthreadingMask = createMask(info->nCorePerComputeUnit, &(info->HyperthreadingMaskWidth));
				coreplusSMT_MaskWidth -= info->HyperthreadingMaskWidth;
      } else {
				info->HyperthreadingMask = 0;
      }
      info->SocketSelectMaskShift = coreplusSMT_MaskWidth;
      info->CoreSelectMask = ( 1 << coreplusSMT_MaskWidth ) - 1; //createMask(coreplusSMT_MaskWidth, &(info->SocketSelectMaskShift));
      info->SocketSelectMaskShift += info->HyperthreadingMaskWidth;
      info->CoreSelectMask <<= info->HyperthreadingMaskWidth;
      info->SocketSelectMask = (-1) ^ (info->CoreSelectMask | info->HyperthreadingMask);

      if ( sched_setaffinity(0, sizeof(cpu_set_t), &BackupMask)) {
				fputs("Error setting affinity!\n",stderr);
				exit(1);
      }
#endif
#ifdef __SYSTEM_DEBUG__
      printf("HTT Mask (%2u)    : ",info->HyperthreadingMaskWidth);
      for (unsigned int j=0; j<32; ++j) {
				if ( info->HyperthreadingMask & (1 << (31 - j)) )
					fputs("1",stdout);
				else
					fputs("0", stdout);
      }
      fputs("\nCore Mask        : ",stdout);
      for (unsigned int j=0; j<32; ++j) {
				if ( info->CoreSelectMask & (1 << (31 - j)) )
					fputs("1",stdout);
				else
					fputs("0", stdout);
      }
      printf("\nSocket Mask (%2u) : ", info->SocketSelectMaskShift);
      for (unsigned int j=0; j<32; ++j) {
				if ( info->SocketSelectMask & (1 << (31 - j)) )
					fputs("1",stdout);
				else
					fputs("0", stdout);
      }
      fputs("\n",stdout);
#endif
      register unsigned int nSockets = 0;
      register unsigned int nCores = 0;
      for (unsigned int thread =0; thread<info->nOverallCores; ++thread) {
				register unsigned int TMask = 1 << (info->APICID[thread] & info->HyperthreadingMask);
				register unsigned int tmp = (info->APICID[thread] & info->CoreSelectMask) >> info->HyperthreadingMaskWidth;
				if (tmp > nCores) nCores = tmp;
				TMask |= 1 << (2 + tmp);

				tmp = (info->APICID[thread] & info->SocketSelectMask) >> info->SocketSelectMaskShift;
				if (tmp > nSockets) nSockets = tmp;
				TMask |= 1 << ( 18 + tmp);
				info->TopologyMasks[thread] = TMask;
#ifdef __SYSTEM_DEBUG__
				fprintf(stderr,"APIC DATA (%2u)   : ", thread);
				for (unsigned int j=0; j<32; ++j) {
					if ( info->APICID[thread] & (1 << (31 - j)) )
						fputs("1",stderr);
					else
						fputs("0", stderr);
				}
				fputs("\n",stderr);
#endif
      }
      //info->nCores = 1 + nCores;
      info->nSockets = 1 + nSockets;
    }
		else {
			fprintf(stderr, "Unrecognized cpu vendor %s\n", info->CPU_Vendor);
			exit(1);
		}
	}

#ifdef __NUMA__
  info->NumaAble     = (numa_available() < 0) ? false : true;
  if (info->NumaAble) {
    info->nNodes       = (unsigned int) numa_num_configured_nodes();
    info->nCpusPerNode = numa_num_configured_cpus() / info->nNodes;
  }
#endif
}

extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
printSystemInfo(const SystemInfo * const info)
{
  char Buffer[80] __attribute__((aligned(16)));
  char OtherBuffer[80] __attribute__((aligned(16)));
  const char Template[80] __attribute__((aligned(16))) = "|                                                                            |\n";
  char * ptr;
  const size_t Length = strlen(Template);
  const char * LastCharacter = &Buffer[Length - 2];

  fputs("|  System informations                                                       |\n"
	"%----------------------------------------------------------------------------%\n", stderr);

#define WRITE(LeftSpace, Title, Value) {\
  memcpy(Buffer, Template, Length*sizeof(char));\
  char * lptr = &Buffer[LeftSpace];\
  lptr += sprintf(lptr, Title);\
  const char * ValuePtr = Value;\
  while (*ValuePtr != '\0' && lptr < LastCharacter) { *lptr++ = *ValuePtr++;}\
  fputs(Buffer,stderr);\
}
  WRITE(5,"Host name      : ", info->Nodename);
  WRITE(5,"User name      : ", info->Username);
  WRITE(5,"Linux kernel   : ", info->Release);
  WRITE(5,"Architecture   : ", info->Architecture);
  WRITE(5,"CPU vendor     : ", info->CPU_Vendor);
  WRITE(5,"CPU Brand      : ", info->CPU_Name);
  sprintf(OtherBuffer, "0x%xh", info->Family);
  WRITE(5,"Family         : ", OtherBuffer);
  sprintf(OtherBuffer, "0x%xh", info->Model);
  WRITE(5,"Model          : ", OtherBuffer);

  ptr = OtherBuffer;
  if (info->Extensions & MM_MMXEXT)   ptr += sprintf(ptr," MMXExt");
  if (info->Extensions & MM_SSE)      ptr += sprintf(ptr," SSE");
  if (info->Extensions & MM_SSE2)     ptr += sprintf(ptr," SSE2");
  if (info->Extensions & MM_SSE3)     ptr += sprintf(ptr," SSE3");
  if (info->Extensions & MM_SSSE3)    ptr += sprintf(ptr," SSSE3");
  if (info->Extensions & MM_SSE41)    ptr += sprintf(ptr," SSE4.1");
  if (info->Extensions & MM_SSE42)    ptr += sprintf(ptr," SSE4.2");
  if (info->Extensions & MM_3DNOW)    ptr += sprintf(ptr," 3DNow");
  if (info->Extensions & MM_3DNOWEXT) ptr += sprintf(ptr," 3DNowExt");
  if (info->Extensions & MM_SSE4A)    ptr += sprintf(ptr," SSE4a");
  if (info->Extensions & MM_POPCOUNT) ptr += sprintf(ptr," POPCOUNT");
  if (info->Extensions & MM_AVX)      ptr += sprintf(ptr," AVX");
  if (ptr > &OtherBuffer[Length-2-21]) {
    ptr = &OtherBuffer[Length-2-21];
    while ( ptr > &OtherBuffer[0] && *ptr != ' ') --ptr;
    *ptr = '\0';
    WRITE(5, "CPU extensions :", OtherBuffer );
    WRITE(5, "               : ", ++ptr);
  } else {
    WRITE(5, "CPU extensions :", OtherBuffer);
  }
  if (info->HyperthreadingAvailable) {
    if ( info->HyperthreadingOn) {
      WRITE(5, "Hyperthreading : ", "Available");
    } else {
      WRITE(5,"Hyperthreading : ", "Available but BIOS disabled");
    }
  }
#if defined(__USE_AFFINITY__) && !defined(_WIN32)
  sprintf(OtherBuffer,"%u", info->nSockets);
  WRITE(5, "Socket         : ", OtherBuffer);
  if (info->nComputeUnit) {
    sprintf(OtherBuffer,"%u", info->nComputeUnit);
    WRITE(5, "Compute unit   : ", OtherBuffer);
    sprintf(OtherBuffer,"%u", info->nCorePerComputeUnit);
    WRITE(5, "Core per unit  : ",OtherBuffer);
  }
  sprintf(OtherBuffer, "%u", info->nCores);
  WRITE(5, "Cores          : ", OtherBuffer);
#endif
  sprintf(OtherBuffer, "%u", info->nOverallCores);
  WRITE(5, "Overall cores  : ", OtherBuffer);
#ifdef __NUMA__
  if (info->NumaAble) {
    WRITE(5,"NUMA           : ", "Available");
    sprintf(OtherBuffer, "%u", info-<nNodes);
    WRITE(5,"NUMA nodes     : ", OtherBuffer);
  } else {
    WRITE(5,"NUMA           : ", "Not available");
  }
#endif

  fputs("%----------------------------------------------------------------------------%\n", stderr);
}

#ifdef _SYSTEM_TEST
int main (int argc, char * argv[])
{
  SystemInfo Info;
  cpu_set_t * Masks;

  getSystemInfo(&Info);
  printSystemInfo(&Info);
#if !defined(_WIN32)
  if (argc == 4) {
    unsigned int SocketId = (unsigned int) atoi(argv[1]);
    unsigned int CoreId = (unsigned int) atoi(argv[2]);
    unsigned int ThreadId = (unsigned int) atoi(argv[3]);
    fprintf(stderr,"Input %u %u %u\n", SocketId, CoreId, ThreadId);

    const unsigned int count = getMasks(&Info, SocketId, CoreId, ThreadId, &Masks);
    printf("%u masks satisfy the criteria\n", count);

    for (unsigned int i=0; i<count; ++i) {
      for (unsigned int j=0; j<32; ++j) {
	if ( CPU_ISSET(31-j,&Masks[i]) )
	  fputs("1",stdout);
	else
	  fputs("0", stdout);
      }
      fputs("\n",stdout);

    }

  }
#endif
  freeSystemInfo(&Info);

  return 1;
}

#endif
#endif
