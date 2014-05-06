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
#include <sys/times.h>
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
    int queryLen = 0;
	int dbLen = 0;
    char query[200];
    char db[200];
    Options options = {-1, -1, DBL_MAX};
    SBMatrix matrix;
    unsigned max_error = atoi( argv[1]);
	
    matrix   = swps3_readSBMatrix( NULL );

    char** read_strs;
    char** ref_strs;
	char* valid_buff;

	read_strs = (char**) malloc(BATCH_RUN * sizeof(char*) );
	ref_strs = (char**) malloc(BATCH_RUN * sizeof(char*) );
	valid_buff = (char*) malloc(BATCH_RUN * sizeof(char) );
	
	int i;
	for (i = 0; i < BATCH_RUN; i++) {
		read_strs[i] = (char*) malloc(200 * sizeof(char) );
		ref_strs[i] = (char*) malloc(200 * sizeof(char) );
	}

    if (argc != 2) {
        printf("$>bin error\n");
        exit(1);
    }

    long long read_size;
    long long read_idx;

	struct tms start_time;
	struct tms end_time;
	struct tms elp_time;

	elp_time.tms_stime = 0;
	elp_time.tms_utime = 0;
	elp_time.tms_cstime = 0;
	elp_time.tms_cutime = 0;

    long long unsigned totalNum = 0;
	long long unsigned passNum = 0;

    int stop = 0;
	do {
        stop = 0;
        for (read_size = 0; read_size < BATCH_RUN; read_size++) {
            fgets(read_strs[read_size], 200, stdin);
            if (strcmp(read_strs[read_size], "end_of_file\n") == 0) {
                stop = 1;
                break;
            }
			int temp_length = strlen(read_strs[read_size]) - 1;
			if (temp_length > queryLen) {
				queryLen = temp_length;
				dbLen = temp_length;
				/* printf("%d\n", queryLen); */
			}

            fgets(ref_strs[read_size], 200, stdin);

			valid_buff[read_size] = 0;
        }

        times(&start_time);

        for (read_idx = 0; read_idx < read_size; read_idx++) {
            memcpy(query, read_strs[read_idx], queryLen);
            swps3_translateSequence(query, queryLen, NULL);
            unsigned score = 0;
            ProfileByte  * profileByte = swps3_createProfileByteSSE( query, queryLen, matrix );
            memcpy(db, ref_strs[read_idx], dbLen);
            swps3_translateSequence(db, dbLen, NULL);

            score = queryLen - swps3_alignmentByteSSE( profileByte, db, dbLen, &options); 
            
			if(score <= max_error)
				valid_buff[read_idx] = 1;
			
			swps3_freeProfileByteSSE( profileByte );
		} 
		
		times(&end_time);

        for (read_idx = 0; read_idx < read_size; read_idx++) {
			if (valid_buff[read_idx]) { 
                passNum++;
                fprintf(stderr, "%s", read_strs[read_idx]);
                fprintf(stderr, "%s", ref_strs[read_idx]);
            }

            totalNum++;
        }

		elp_time.tms_stime += end_time.tms_stime - start_time.tms_stime;
		elp_time.tms_utime += end_time.tms_utime - start_time.tms_utime;
		
		if (stop)
			break;

    } while (!stop);
	
	fprintf(stderr, "end_of_file\n");
	printf("passNum:\t%lld\n", passNum);
	printf("totalNum:\t%lld\n", totalNum);
	printf("total_time: %f\n", (double) elp_time.tms_utime / sysconf(_SC_CLK_TCK) ); 

	for (i = 0; i < BATCH_RUN; i++) {
		free(read_strs[i]);
		free(ref_strs[i]);
	}
	free(read_strs);
	free(ref_strs);
	free(valid_buff);


    return 0;
}

