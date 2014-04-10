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
#include <sys/time.h>
#include <string.h>
#define TEST_SIZE 1
#define LOC_SIZE 1

#define BATCH_RUN 1000000

static ssize_t write_all(int fd, const void *data, size_t len) {
    size_t sent = 0;
    while(sent < len) {
        ssize_t res;
        res = write(fd,(const int8_t*)data+sent,len-sent);
        switch(res) {
            case 0:
                return sent;
            case -1:
                return -1;
            default:
                sent += res;
        }
    }
    return sent;
}

static ssize_t read_all(int fd, void *buf, size_t len) {
    size_t recv = 0;
    while(recv < len) {
        ssize_t res;
        res = read(fd,(int8_t*)buf+recv,len-recv);
        switch(res) {
            case 0:
                return recv;
            case -1:
                return -1;
            default:
                recv += res;
        }
    }
    return recv;
}

int main( int argc, char * argv[] ){
    char * matrixFile = NULL, * queryFile = NULL, * dbFile = NULL;
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

    matrix   = swps3_readSBMatrix( matrixFile );
            int *pipes_read;
            int *pipes_write;
            char **seq_names;
            pid_t *children;
            int child_id = -1;
            int childpipe_read = -1;
            int childpipe_write = -1;
            FastaLib * dbLib = swps3_openLib( dbFile );
#if defined(__SSE2__)
            ProfileByte  * profileByte = swps3_createProfileByteSSE( query, queryLen, matrix );
            ProfileShort * profileShort = swps3_createProfileShortSSE( query, queryLen, matrix );
#endif
#if defined(__PS3)
            SPEProfile * profileByte = swps3_createProfileBytePPU(query, queryLen, matrix, MAX_SEQ_LENGTH);
            SPEProfile * profileShort = swps3_createProfileShortPPU(query, queryLen, matrix, MAX_SEQ_LENGTH);
            SPEProfile * profileFloat = swps3_createProfileFloatPPU(query, queryLen, matrix, MAX_SEQ_LENGTH);
            /* by default byte profile will be loaded */
            int current_profile_is_byte = 0;
            /* loadProfileByte(profileByte, MAX_SEQ_LENGTH, &options);*/
#endif
#if defined(__ALTIVEC__)
            void *profileByteAltivec = swps3_createProfileByteAltivec(query, queryLen, matrix);
            void *profileShortAltivec = swps3_createProfileShortAltivec(query, queryLen, matrix);
            void *profileFloatAltivec = swps3_createProfileFloatAltivec(query, queryLen, matrix);
#endif
            pipes_read = malloc(threads*sizeof(*pipes_read));
            pipes_write = malloc(threads*sizeof(*pipes_write));
            children = malloc(threads*sizeof(*children));
            seq_names = malloc(threads*sizeof(*seq_names));
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


    tms start_time;
    tms end_time;
    tms elp_time;

    elp_time.tms_stime = 0;
    elp_time.tms_utime = 0;
    elp_time.tms_cstime = 0;
    elp_time.tms_cutime = 0;

    int stop = 0;

    while (!stop) {  
        stop = 0;
        for (read_size = 0; read_size < BATCH_RUN; read_size++) {		
            fgets(read_strs[read_size], 200, stdin)
                if (strcmp(read_strs[read_size], "end_of_file\0") == 0) {
                    stop = 1;
                    break;
                }

            fgets(ref_strs[read_size], 200, stdin)
        }

        times(&start_time);

        for (read_idx = 0; read_idx < read_size; read_idx++) {
            /* TODO: load query */
            query = read_strs[read_idx];
            queryLen = strlen(read_strs[read_idx]);
            double score = 0;

                int dbLen;
                char * db;
                char * dbName;

                db = ref_strs[read_idx];
                dbLen = strlen(ref_strs[read_idx]);
                strncpy(dbName,itoa(read_idx),MAX_SEQ_NAME_LENGTH);
                if(db == NULL) break;

#ifdef DEBUG
                for(i=0; i<queryLen; ++i) printf("\t%c",query[i]);
                printf("\n");
#endif

#if defined(__SSE2__)
                if(type == SSE2) {
                    if( (score = swps3_alignmentByteSSE( profileByte, db, dbLen, &options )) >= DBL_MAX ) {
                        score = swps3_alignmentShortSSE( profileShort, db, dbLen, &options );
                        assert(score >= 250 && "score too low");
                    }
                }
#endif
#if defined(__PS3)
                if(type == PS3) {
#if defined(__ALTIVEC__)
                    if(child_id == 6) {
#if 0
                        score = swps3_dynProgrFloatAltivec(db, dbLen, profileFloatAltivec, &options);
#else
                        score = swps3_dynProgrByteAltivec(db, dbLen, profileByteAltivec, &options);
                        if(score >= DBL_MAX)
                            score = swps3_dynProgrShortAltivec(db, dbLen, profileShortAltivec, &options);
#endif
                    } else {
#endif /* __ALTIVEC__ */
#if 0
                        loadProfileFloat(profileFloat, MAX_SEQ_LENGTH, &options);
                        score = alignmentProfileSPE( db, dbLen );
#elif 1
                        if(!current_profile_is_byte) swps3_loadProfileByte(profileByte, MAX_SEQ_LENGTH, &options);
                        current_profile_is_byte = 1;
                        score = swps3_alignmentProfileSPE( db, dbLen );
                        if( score >= DBL_MAX ) {
                            if(current_profile_is_byte) swps3_loadProfileShort(profileShort, MAX_SEQ_LENGTH, &options);
                            current_profile_is_byte = 0;
                            score = swps3_alignmentProfileSPE( db, dbLen );
                        }
#elif 1
                        if(!current_profile_is_byte) swps3_loadProfileShort(profileShort, MAX_SEQ_LENGTH, &options);
                        current_profile_is_byte = 1;
                        score = swps3_alignmentProfileSPE( db, dbLen );
#else
                        score = swps3_alignmentByteSPE( matrix, query, queryLen, db, dbLen, &options );
                        if( score >= DBL_MAX ) {
                            score = swps3_alignmentShortSPE( matrix, query, queryLen, db, dbLen, &options );
                        }
#endif
#if defined(__ALTIVEC__)
                    }
#endif /* __ALTIVEC__ */
                }
#endif /* __PS3 */
#if defined(__ALTIVEC__)
                if(type == ALTIVEC) {
#if 0
                    score = swps3_dynProgrFloatAltivec(db, dbLen, profileFloatAltivec, &options);
#else
                    score = swps3_dynProgrByteAltivec(db, dbLen, profileByteAltivec, &options);
                    if(score >= DBL_MAX)
                        score = swps3_dynProgrShortAltivec(db, dbLen, profileShortAltivec, &options);
                }
#endif
#endif /* __ALTIVEC__ */
                if(type == SCALAR) {
                    double dmatrix[MATRIX_DIM*MATRIX_DIM];
                    for(i=0;i<MATRIX_DIM*MATRIX_DIM;++i) dmatrix[i]=matrix[i];
                    score = swps3_alignScalar( dmatrix, query, queryLen, db, dbLen, &options);
                }
            }

            if(childpipe_write > 0) {
                ssize_t res;

                res = write_all(childpipe_write,&score,sizeof(score));
                if(res != sizeof(score)) {
                    perror("error during write()");
                    exit(1);
                }
            } else {
                if(score >= options.threshold) {
                    /* printf(">threshold\t%s\n",dbName); */
                } else {
                    /* printf("%g\t%s\n",score,dbName); */
                    align_num++;
                    fprintf(stderr, "%s\n", read_strs[read_idx] );
                    fprintf(stderr, "%s\n", ref_strs[read_idx] );
                }
            }

            dCount++; dResidues+=dbLen;

            total_num++;
    }


#if defined(__SSE2__)
    swps3_freeProfileByteSSE( profileByte );
    swps3_freeProfileShortSSE( profileShort );
#endif
#if defined(__PS3)
    swps3_freeProfilePPU( profileByte );
    swps3_freeProfilePPU( profileShort );
    swps3_freeProfilePPU( profileFloat );
#endif
    swps3_closeLib( dbLib );
}
fprintf(stderr,"%d[%d] x %d[%d]\n", qCount, qResidues, dCount, dResidues );

/* swps3_closeLib( queryLib ); */
return 0;
}
