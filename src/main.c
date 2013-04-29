// Copyright (c) 2013 Genome Research Limited.
//
// This file is part of samtools merge.
//
// samtools merge is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// this program. If not, see L<http://www.gnu.org/licenses/>.

#define _GNU_SOURCE

#include <htslib/sam.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>

struct parsed_opts {
    size_t input_count;
    char** input_name;
    char* output_name;
};

struct state {
    size_t input_count;
    samFile** input_file;
    bam_hdr_t** input_header;
    samFile* output_file;
    bam_hdr_t* output_header;
};

typedef struct parsed_opts parsed_opts_t;
typedef struct state state_t;

parsed_opts_t* parse_args(int argc, char** argv) {
    if (argc < 3) {
        dprintf(STDERR_FILENO, "Arguments should be: merge <input1.bam> [<inputX.bam> ...] <output.bam>\r\n");
        return NULL;
    }
    
    parsed_opts_t* retval = malloc(sizeof(parsed_opts_t));
    if (! retval ) return NULL;

    retval->input_count = argc-2;
    retval->input_name = (char**)calloc(retval->input_count,sizeof(char*));
    size_t i = 0;
    for (; i < retval->input_count; i++) {
        retval->input_name[i] = strdup(argv[i+1]);
    }

    retval->output_name = strdup(argv[i+1]);

    return retval;
}

bam_hdr_t* merge_headers( const bam_hdr_t** input_header, const size_t input_count ) {
    if (input_count == 0) return NULL;
    bam_hdr_t* retval = (bam_hdr_t*)malloc(sizeof(bam_hdr_t));
    // TODO: merge headers instead of just taking first one
    // TODO: Need to clone interior of this?
    memcpy((void*)retval, (const void*)input_header[0], sizeof(bam_hdr_t));
    return retval;
}

bool init(parsed_opts_t* opts, state_t** state_out) {
    state_t* retval = (state_t*) malloc(sizeof(state_t));
    if (retval == NULL) {
        dprintf(STDERR_FILENO, "Out of memory");
        return false;
    }
    *state_out = retval;
    retval->input_count = opts->input_count;
    
    // Open files
    retval->input_file = (samFile**)calloc(opts->input_count, sizeof(samFile*));
    retval->input_header = (bam_hdr_t**)calloc(opts->input_count, sizeof(bam_hdr_t*));
    for (size_t i = 0; i < opts->input_count; i++) {
        retval->input_file[i] = sam_open(opts->input_name[i], "rb", 0);
        if (retval->input_file[i] == NULL) {
            dprintf(STDERR_FILENO, "Could not open input file: %s\r\n", opts->input_name[i]);
            return false;
        }
        retval->input_header[i] = sam_hdr_read(retval->input_file[i]);
    }
    
    retval->output_header = merge_headers(retval->input_header, opts->input_count);
    // TODO: create SQ translation table
    // TODO: create RG translation table

    retval->output_file = sam_open(opts->output_name, "wb", 0);
    
    if (retval->output_file == NULL) {
        dprintf(STDERR_FILENO, "Could not open output file: %s\r\n", opts->output_name);
        return false;
    }
    
    return true;
}

size_t selectRead( bam1_t **file_read, size_t input_count )
{
    assert(input_count != 0);
    // need to find element with lowest tid and pos
    size_t min = SIZE_T_MAX;
    // horrible hack of the day
    // treat the int32_t tid as a uint32_t so that -1 (aka unmapped) is treated as UINT32_MAX
    uint32_t tid_min;
    int32_t pos_min;
    // load initial value
    size_t i = 0;
    for (; i < input_count; i++) {
        if (file_read[i] != NULL) {
            tid_min = (uint32_t)file_read[i]->core.tid;
            pos_min = file_read[i]->core.pos;
            min = i;
            break;
        }
    }
    assert(min != SIZE_T_MAX); // No valid files?
    
    // then resume our search
    for (;i < input_count; i++) {
        if (file_read[i] != NULL) {
            // To complicate matters tid == -1 is a special value which should always go last
            if ((tid_min > (uint32_t)file_read[i]->core.tid ) ||
                (tid_min == (uint32_t)file_read[i]->core.tid && pos_min > file_read[i]->core.pos)) {
                tid_min = (uint32_t)file_read[i]->core.tid;
                pos_min = file_read[i]->core.pos;
                min = i;
            }
        }
    }
    
    assert(min != SIZE_T_MAX); // No valid files?
    
    return min;
}

bool merge(state_t* state) {
    if (sam_hdr_write(state->output_file, state->output_header) != 0) {
        return false;
    }
    
    bam1_t** file_read = calloc(state->input_count, sizeof(bam1_t*));
    size_t files_to_merge = state->input_count;
    for (size_t i = 0; i < state->input_count; i++) {
        file_read[i] = bam_init1();
        if (sam_read1(state->input_file[i], state->input_header[i], file_read[i]) < 0) {
            bam_destroy1(file_read[i]);
            file_read[i] = NULL;
            files_to_merge--;
        }
    }

    while (files_to_merge > 0) {
        size_t i = selectRead(file_read, state->input_count);
        sam_write1(state->output_file, state->output_header, file_read[i]);
        if (sam_read1(state->input_file[i], state->input_header[i], file_read[i]) < 0) {
            bam_destroy1(file_read[i]);
            file_read[i] = NULL;
            files_to_merge--;
        }
    }

    // Clean up
    for (size_t i = 0; i < state->input_count; i++) { if (file_read[i]) { bam_destroy1(file_read[i]); } }

    return true;
}

void cleanup_opts(parsed_opts_t* opts) {
    free(opts->output_name);
    for (size_t i = 0; i < opts->input_count; i++) {
        free(opts->input_name[i]);
    }
}

void cleanup_state(state_t* state) {
    sam_close(state->output_file);
    bam_hdr_destroy(state->output_header);
    for (size_t i = 0; i < state->input_count; i++) {
        sam_close(state->input_file[i]);
        bam_hdr_destroy(state->input_header[i]);
    }
}

int main(int argc, char** argv) {

    parsed_opts_t* opts = parse_args(argc, argv);
    state_t* state = NULL;
    if (!opts || !init(opts, &state)) return -1;
    
    if (!merge(state)) return -1;
    
    cleanup_opts(opts);
    cleanup_state(state);
    
    return 0;
}
