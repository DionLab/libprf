#include "prf_config.h"
#include "java/PfTools_Profile.h"
#include "java/PfTools_Matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#define __USE_INLINE_FUNCTIONS__
#include "pfProfile.h"
#include "pfSequence.h"
#include "pfCompute.h"
#include "pfOutput.h"

/*
 * Class:     Profile
 * Method:    loadProfile
 * Signature: (Ljava/lang/String;)J
 */
JNIEXPORT jboolean JNICALL Java_PfTools_Profile_loadProfile(JNIEnv *jEnv, jobject ThisObj,
                                                            jstring ProfileFileName)
{
	/* Allocate memory for structure */
	struct Profile * prf = malloc(sizeof(struct Profile));
	if (prf == NULL) {
		fputs("Memory allocation of profile structure failed\n", stderr);
		return JNI_FALSE ;
	}
	/* Retrieve File name */
	const char *prfName = (*jEnv)->GetStringUTFChars(jEnv, ProfileFileName, NULL);
	if (NULL == prfName) {
		fputs("Unable to get profile file name from java\n", stderr);
		return JNI_FALSE;
	}
	printf("Profile file name : %s\n", prfName);
	
	/* Read profile */
	const int ReadErr = ReadProfile(prfName, prf, true, false);
	(*jEnv)->ReleaseStringUTFChars(jEnv, ProfileFileName, prfName); 
	if (ReadErr != 1) {
		fprintf(stderr, "Profile reading returned %i!\n", ReadErr);
		return JNI_FALSE; 
	}
	
	/* Reserve memory for work array */
	const size_t WorkSize = 4UL*(prf->Length + 1UL)*sizeof(int) + 63UL;
	int * WORK = _mm_malloc(WorkSize, 64);
	if (WORK == NULL) {
		FreeProfile(prf, true);
		return JNI_FALSE;
	}
	
	/* Set profile level to 0, mode 1 */
	if (SetProfileLevelAndMode(prf, 0, 1) != 0) {
		fprintf(stderr, "Unable to get profile Level 0, mode 1\n");
		return JNI_FALSE;
	}
	
	/*Set data to Java */
  jclass thisClass = (*jEnv)->GetObjectClass(jEnv, ThisObj);
  jfieldID fidNumber = (*jEnv)->GetFieldID(jEnv, thisClass, "profilePointer", "J");
  if (NULL == fidNumber) {
		fputs("Unable to get field id for profile pointer\n", stderr);
		goto bail;
	}
  (*jEnv)->SetLongField(jEnv, ThisObj, fidNumber, (jlong) (uintptr_t) prf);
	fidNumber = (*jEnv)->GetFieldID(jEnv, thisClass, "workPointer", "J");
  if (NULL == fidNumber) {
		fputs("Unable to get field id for work pointer\n", stderr);
		goto bail;
	}
  (*jEnv)->SetLongField(jEnv, ThisObj, fidNumber, (jlong) (uintptr_t) WORK);
	
	return JNI_TRUE;
	
	bail:;
	_mm_free(WORK);
	FreeProfile(prf, true);
	return JNI_FALSE ;
}

/*
 * Class:     Profile
 * Method:    freeProfile
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_PfTools_Profile_freeProfile(JNIEnv *jEnv, jobject ThisObj)
{
	/* Get back profile and Work */
	jclass thisClass = (*jEnv)->GetObjectClass(jEnv, ThisObj);
  jfieldID fidNumber = (*jEnv)->GetFieldID(jEnv, thisClass, "profilePointer", "J");
  if (NULL != fidNumber) {
		struct Profile * const prf = (struct Profile*) (*jEnv)->GetLongField(jEnv, ThisObj, fidNumber);
		if (prf) FreeProfile(prf, true);
	}
	fidNumber = (*jEnv)->GetFieldID(jEnv, thisClass, "workPointer", "J");
  if (NULL != fidNumber) {
		int * const WORK = (int*) (*jEnv)->GetLongField(jEnv, ThisObj, fidNumber);
		if (WORK) _mm_free(WORK);
	}
}

/*
 * Class:     Profile
 * Method:    createJSON
 * Signature: (Ljava/lang/String;)V
 */
JNIEXPORT jboolean JNICALL Java_PfTools_Profile_createJSON(JNIEnv *jEnv, jobject ThisObj, jstring FileName)
{
	jboolean ret = JNI_FALSE;
	
	/* Get back profile  */
	jclass thisClass = (*jEnv)->GetObjectClass(jEnv, ThisObj);
  jfieldID fidNumber = (*jEnv)->GetFieldID(jEnv, thisClass, "profilePointer", "J");
  if (NULL == fidNumber) goto bail;
  const struct Profile * const prf = (struct Profile*) (*jEnv)->GetLongField(jEnv, ThisObj, fidNumber);
	
	/* Retrieve File name */
	const char *FName = (*jEnv)->GetStringUTFChars(jEnv, FileName, NULL);
	if (NULL == FName) {
		fputs("Unable to get JSON file name from java\n", stderr);
		goto bail;
	}
	
	FILE *out = fopen(FName, "w");
	if (out) {
		PrintJSON(prf, out);
		fclose(out);
		ret = JNI_TRUE;
	}
		
	(*jEnv)->ReleaseStringUTFChars(jEnv, FileName, FName);
	bail:;
		return ret;
}

/*
 * Class:     Matrix
 * Method:    computeMatrix
 * Signature: (Ljava/lang/String;)V
 */
JNIEXPORT jboolean JNICALL Java_PfTools_Matrix_buildMatrix(JNIEnv *jEnv, jobject ThisObj, jstring Sequence)
{
	union lScores * Matrix = NULL;
	/* Retrieve Sequence */
	const char *Seq = (*jEnv)->GetStringUTFChars(jEnv, Sequence, NULL);
	if (NULL == Seq) {
		fputs("Unable to accecc sequence\n", stderr);
		return JNI_FALSE;
	}
	/* Get back profile and Work */
	jclass thisClass = (*jEnv)->GetObjectClass(jEnv, ThisObj);
  jfieldID fidNumber = (*jEnv)->GetFieldID(jEnv, thisClass, "profilePointer", "J");
  if (NULL == fidNumber) {
		fputs("Unable to get field id for profile pointer\n", stderr);
		goto bail;
	}
  const struct Profile * const prf = (struct Profile*) (*jEnv)->GetLongField(jEnv, ThisObj, fidNumber);
	if (prf == NULL) {
		fputs("Invalid profile pointer retrieved\n", stderr);
		goto bail;
	}
//  	struct Profile * prf = malloc(sizeof(struct Profile));
//  	if (ReadProfile("/tmp/junit4696355374887423/sh3.prf", prf, true, false) !=1) {
// 		fputs("NO PROFILE!\n", stderr);
// 		exit(1);
// 	}
	fidNumber = (*jEnv)->GetFieldID(jEnv, thisClass, "workPointer", "J");
  if (NULL == fidNumber) {
		fputs("Unable to get field id for work pointer\n", stderr);
		goto bail;
	}
	int * const WORK = (int*) (*jEnv)->GetLongField(jEnv, ThisObj, fidNumber);
	if (WORK == NULL) {
		fputs("Invalid work pointer retrieved\n", stderr);
		goto bail;
	}
	
	const size_t SeqLength = strlen(Seq);
	
	if (SeqLength < prf->Length) {
		fputs("Sequence has to be at least the size of profile here\n", stderr);
		return JNI_FALSE;
	}
 	const size_t MatrixSize = (SeqLength+1UL)*(prf->Length+1UL)*sizeof(union lScores);
	Matrix = _mm_malloc(MatrixSize, 16);
	if (Matrix == NULL) {
		fputs("Invalid matrix pointer retrieved\n", stderr);
		goto bail;
	}
	
	unsigned char * const SequenceIndex = malloc(SeqLength*sizeof(unsigned char));
	if (SequenceIndex == NULL) {
			fputs("Unable to allocate memory for sequence index\n", stderr);
			_mm_free(Matrix);
			return JNI_FALSE;
	}
	
	PFSequence PFSeq = { .ProfileIndex=SequenceIndex , .Length = SeqLength} ;
	
	
	/* TODO: Remove that memory copy and perform both in a row */
	/* Copy sequence to local space index */ 
	memcpy(PFSeq.ProfileIndex, Seq, PFSeq.Length);
	
	/* Translate into indices */
	TranslateSequenceToIndex(&PFSeq, prf->Alphabet_Mapping);

	
	/* Compute matrix */
	FPGA_sse41.BuildMatrix(prf, PFSeq.ProfileIndex, Matrix, WORK, NULL, 0, SeqLength);
// 	memset(Matrix, 0, MatrixSize);
	
// 	int i=0;
// 	printf("Input sequence of length %lu\n", SeqLength);
// 	while (i < SeqLength) { printf("%.60s\n", &Seq[i]); i += 60; }
// 	
// 	struct IO_Data Options = { .IsBinary = false, .Separate = false };
// 	Options.BaseFileName = "/tmp/oops";
// 	fprintf(stderr, "About to dump to %s\n", Options.BaseFileName);
// 	TabOutput(Matrix, NULL, Seq, "Test", prf, &Standard_sse41, SeqLength, &Options, NULL);
// 	
	
	/* Free sequence index */
	free(SequenceIndex);
	
	/* Set data to Java */
	fidNumber = (*jEnv)->GetFieldID(jEnv, thisClass, "matrixPointer", "J");
  if (NULL == fidNumber) {
		fputs("Unable to get field id for matrix pointer\n", stderr);
		goto bail;
	}
	(*jEnv)->SetLongField(jEnv, ThisObj, fidNumber, (jlong) (uintptr_t) Matrix);
	fidNumber = (*jEnv)->GetFieldID(jEnv, thisClass, "SequenceLength", "J");
  if (NULL == fidNumber) {
		fputs("Unable to get field id for sequence length\n", stderr);
		goto bail;
	}
	(*jEnv)->SetLongField(jEnv, ThisObj, fidNumber, (jlong) SeqLength);
	
	(*jEnv)->ReleaseStringUTFChars(jEnv, Sequence, Seq);
	return JNI_TRUE;
	
	bail:;
		(*jEnv)->ReleaseStringUTFChars(jEnv, Sequence, Seq);
		if (Matrix) _mm_free(Matrix);
		return JNI_FALSE;
}

/*
 * Class:     Matrix
 * Method:    dumpTabulated
 * Signature: (Ljava/lang/String;)V
 */
JNIEXPORT void JNICALL Java_PfTools_Matrix_dumpTabulated(JNIEnv *jEnv, jobject ThisObj, jstring FileName)
{
	const static char Header[] = "JNI_Matrix_output"; 
	struct IO_Data Options = { .IsBinary = false, .Separate = false };
	
	/* Get back profile and Work */
	jclass thisClass = (*jEnv)->GetObjectClass(jEnv, ThisObj);
  jfieldID fidNumber = (*jEnv)->GetFieldID(jEnv, thisClass, "profilePointer", "J");
  if (NULL == fidNumber)  {
		fputs("Error getting profile structure member\n", stderr);
		goto bail;
	}
  const struct Profile * const prf = (struct Profile*) (*jEnv)->GetLongField(jEnv, ThisObj, fidNumber);
	if (prf == NULL) {
		fputs("Error getting profile structure pointer\n", stderr);
		goto bail;
	}
	fidNumber = (*jEnv)->GetFieldID(jEnv, thisClass, "matrixPointer", "J");
  if (NULL == fidNumber) {
		fputs("Error getting matrix array member\n", stderr);
		goto bail;
	}
	union lScores * const Matrix = (union lScores*) (*jEnv)->GetLongField(jEnv, ThisObj, fidNumber);
	if (NULL == Matrix) {
		fputs("Error getting matrix array pointer\n", stderr);
		goto bail;
	}
	fidNumber = (*jEnv)->GetFieldID(jEnv, thisClass, "SequenceLength", "J");
  if (NULL == fidNumber)  {
		fputs("Error getting sequence length member\n", stderr);
		goto bail;
	}
	const size_t SeqLength =  (size_t) (*jEnv)->GetLongField(jEnv, ThisObj, fidNumber);
	
	Options.BaseFileName = (char*) (*jEnv)->GetStringUTFChars(jEnv, FileName, NULL);
	if (Options.BaseFileName == NULL) {
		fputs("Error reding tabulated output file name\n", stderr);
		goto bail;
	}
	
	if (SeqLength) {
		fprintf(stderr, "About to run Tab output to base name %s\n", Options.BaseFileName);
		TabOutput(Matrix, NULL, NULL, Header, prf, &FPGA_sse41, SeqLength, &Options, NULL);
	}
	else {
			fputs("Sequence length is null!\n", stderr);
	}
	(*jEnv)->ReleaseStringUTFChars(jEnv, FileName, Options.BaseFileName);
	
	bail:;
		return;
}

/*
 * Class:     Matrix
 * Method:    freeMatrix
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_PfTools_Matrix_freeMatrix(JNIEnv *jEnv, jobject ThisObj)
{
	/* Get back the matrix */
	jclass thisClass = (*jEnv)->GetObjectClass(jEnv, ThisObj);
	jfieldID fidNumber = (*jEnv)->GetFieldID(jEnv, thisClass, "matrixPointer", "J");
	if (NULL != fidNumber) {
		union lScores * Matrix = (union lScores*)  (*jEnv)->GetLongField(jEnv, ThisObj, fidNumber);
		if (Matrix) {
			_mm_free(Matrix);
			(*jEnv)->SetLongField(jEnv, ThisObj, fidNumber, (jlong) 0);
		}
	}
}

