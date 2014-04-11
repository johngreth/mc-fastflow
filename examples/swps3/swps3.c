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
#include "DynProgr_sse_byte.h"
#include "DynProgr_sse_short.h"
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

#define BATCH_RUN 10000

int main( int argc, char * argv[] ){
    printf("start\n");
    fflush(stdout);
    int queryLen;
    char * query;
    Options options = {-12,-2,DBL_MAX};
    int qCount=0, dCount=0, qResidues=0, dResidues=0;
    SBMatrix matrix;
    options.gapOpen = 1;
    options.gapExt = 1;
    options.threshold = atoi( argv[1]);
    printf("Pre matrix\n");
    matrix   = swps3_readSBMatrix( NULL );
    printf("Post matrix\n");
    char read_strs[BATCH_RUN][200];
    char ref_strs[BATCH_RUN][200];

    if (argc != 2) {
        printf("$>bin error\n");
        exit(1);
    }
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
            if (strcmp(read_strs[read_size], "end_of_file\n") == 0) {
                stop = 1;
                break;
            }

            fgets(ref_strs[read_size], 200, stdin);
        }

        start = clock();

        for (read_idx = 0; read_idx < read_size; read_idx++) {
            query = read_strs[read_idx];
            queryLen = strlen(read_strs[read_idx]);
            double score = 0;
            ProfileByte  * profileByte = swps3_createProfileByteSSE( query, queryLen, matrix );
            ProfileShort * profileShort = swps3_createProfileShortSSE( query, queryLen, matrix );
            qCount++; qResidues+=queryLen;
            dCount=dResidues=0;

            int dbLen;
            char * db;

            db = ref_strs[read_idx];
            dbLen = strlen(ref_strs[read_idx]);
            if(db == NULL) break;

                if( (score = swps3_alignmentByteSSE( profileByte, db, dbLen, &options )) >= DBL_MAX ) {
                    score = swps3_alignmentShortSSE( profileShort, db, dbLen, &options );
                    /* assert(score >= 250 && "score too low"); */
                }
            if(score < options.threshold) {
                align_num++;
                fprintf(stderr, "%s", read_strs[read_idx] );
                fprintf(stderr, "%s", ref_strs[read_idx] );
            }

            dCount++; dResidues+=dbLen;

            total_num++;
        }
        total += clock() - start;
    }
    fprintf(stderr,"%d[%d] x %d[%d]\n", qCount, qResidues, dCount, dResidues );

    return 0;
}

