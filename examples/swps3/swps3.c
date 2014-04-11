/** \file swps3.c
 *
 * Main procedure and multi-threading code.
 */
/*
 * Copyright (c) 2007-2008 ETH Zürich, Institute of Computational Science
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "swps3.h"
#include "matrix.h"
#include "fasta.h"
#include "DynProgr_scalar.h"
#if defined(__SSE2__)
#include "DynProgr_sse_byte.h"
#include "DynProgr_sse_short.h"
#endif
#if defined(__ALTIVEC__)
#include "DynProgr_altivec.h"
#endif
#if defined(__PS3)
#include "DynProgr_PPU.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/select.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <time.h>
#include <string.h>
#define TEST_SIZE 1
#define LOC_SIZE 1

#define BATCH_RUN 1000000

int main( int argc, char * argv[] ){
    printf("start\n");
    fflush(stdout);
    char * matrixFile = NULL;
    int i, queryLen;
    char * query;
#if defined(__SSE2__)
    SWType type = SSE2;
#elif defined(__PS3)
    SWType type = PS3;
#elif defined(__ALTIVEC__)
    SWType type = ALTIVEC;
#else
    SWType type = SCALAR;
#endif
    Options options = {-12,-2,DBL_MAX};
    int qCount=0, dCount=0, qResidues=0, dResidues=0;
#ifdef HAVE_SYSCONF_NPROCESSORS
    int threads = sysconf(_SC_NPROCESSORS_ONLN);
#else
#if defined(__PS3)
    int threads = 6;
#else
    int threads = 1;
#endif
#endif
    SBMatrix matrix;
    /*        FastaLib * queryLib;
              for ( i=1; i<argc; i++ ){
              if (argv[i][0]=='-'){
              switch( argv[i][1] ){
              case 'h':
              matrixFile = NULL;
              i = argc; break;
              case 's':
              type = SCALAR;
              break;
              case 't':
              options.threshold = atoi( argv[++i] );
              break;
              case 'i':
              options.gapOpen = atoi( argv[++i] );
              break;
              case 'e':
              options.gapExt = atoi( argv[++i] );
              break;
              case 'j':
              threads = atoi( argv[++i] );
              break;
              default:
              matrixFile = NULL;
              i = argc; break;
              }
              }else{
              if (matrixFile == NULL)
              matrixFile = argv[i];
              else if (queryFile == NULL)
              queryFile = argv[i];
              else if (dbFile == NULL)
              dbFile = argv[i];
              else{
              matrixFile = NULL;
              i = argc; break;
              }
              }
              }
              if ( matrixFile == NULL || queryFile == NULL || dbFile == NULL ){
              printf( "Usage: %s [-h] [-s] [-j num] [-i num] [-e num] [-t num] matrix query db\n", argv[0] );
              return 0;
              }
              */
    options.gapOpen = 1;
    options.gapExt = 1;
    threads = 1;
    options.threshold = atoi( argv[1]);
    printf("Pre matrix\n");
    matrix   = swps3_readSBMatrix( matrixFile );
    printf("Post matrix\n");
    /* queryLib = swps3_openLib( queryFile );

       START BENCHMARK CODE
       */
    char read_strs[BATCH_RUN][200];
    char ref_strs[BATCH_RUN][200];

    if (argc != 2) {
        printf("$>bin error\n");
        exit(1);
    }

    /* FILE *input;
       FILE *output;

       input = fopen(argv[1], "r");
       output = fopen(argv[2], "w");
       */

    long long align_num = 0;
    long long total_num = 0;

    long long read_size;
    long long read_idx;


    clock_t start;
    long long unsigned total;

    int stop = 0;
printf("loop\n");
    while (!stop) {  
        stop = 0;
        for (read_size = 0; read_size < BATCH_RUN; read_size++) {		
            fgets(read_strs[read_size], 200, stdin);
            if (strcmp(read_strs[read_size], "end_of_file\0") == 0) {
                stop = 1;
                break;
            }

            fgets(ref_strs[read_size], 200, stdin);
        }

        start = clock();

        for (read_idx = 0; read_idx < read_size; read_idx++) {
            /* TODO: load query */
            query = read_strs[read_idx];
            queryLen = strlen(read_strs[read_idx]);
            double score = 0;

/*
    int *pipes_read;
    int *pipes_write;
    char **seq_names;
    pid_t *children;
    int child_id = -1;
    int childpipe_read = -1;
    int childpipe_write = -1;
*/

#if defined(__SSE2__)
    ProfileByte  * profileByte = swps3_createProfileByteSSE( query, queryLen, matrix );
    ProfileShort * profileShort = swps3_createProfileShortSSE( query, queryLen, matrix );
#endif
/*
    pipes_read = malloc(threads*sizeof(*pipes_read));
    pipes_write = malloc(threads*sizeof(*pipes_write));
    children = malloc(threads*sizeof(*children));
    seq_names = malloc(threads*sizeof(*seq_names));
*/
    for(i=0;i<threads;++i) {
        pipes_read[i]=-1;
        pipes_write[i]=-1;
        children[i]=-1;
        seq_names[i]=malloc((MAX_SEQ_NAME_LENGTH+1)*sizeof(char));
        seq_names[i][MAX_SEQ_NAME_LENGTH]='\0';
    }

    qCount++; qResidues+=queryLen;
    dCount=dResidues=0;

    childpipe_read = -1;
    childpipe_write = -1;
    child_id = 0;
            int dbLen;
            char * db;
            char dbName[30];

            db = ref_strs[read_idx];
            dbLen = strlen(ref_strs[read_idx]);
            sprintf(dbName, "%d", (int)read_idx);
            if(db == NULL) break;

#ifdef DEBUG
            for(i=0; i<queryLen; ++i) printf("\t%c",query[i]);
            printf("\n");
#endif

            if(type == SSE2) {
                if( (score = swps3_alignmentByteSSE( profileByte, db, dbLen, &options )) >= DBL_MAX ) {
                    score = swps3_alignmentShortSSE( profileShort, db, dbLen, &options );
                    assert(score >= 250 && "score too low");
                }
            }
            if(score < options.threshold) {
                align_num++;
                fprintf(stderr, "%s\n", read_strs[read_idx] );
                fprintf(stderr, "%s\n", ref_strs[read_idx] );
            }

            total += clock() - start;

            dCount++; dResidues+=dbLen;

            total_num++;
        }
    }
    fprintf(stderr,"%d[%d] x %d[%d]\n", qCount, qResidues, dCount, dResidues );

    return 0;
}

