/*******************************************************
                        PFTOOLS
 *******************************************************
  Aug 23, 2013 version.c
 *******************************************************
 (C) 2013 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/

#include "prf_config.h"
#include "pfVersion.h"

char Version[] = PRF_VERSION;

const char * GetVersion() 
{
  return Version;
}
