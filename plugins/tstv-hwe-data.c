/*  plugins/counts.c -- counts SNPs, Indels, and total number of sites.

    Copyright (C) 2013, 2014 Genome Research Ltd.

Author: Petr Danecek <pd3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>

typedef struct
{
  int argc;
  char **argv;
  bcf_hdr_t *hdr;
  int32_t *gt_arr;
  int mgt_arr;
}
args_t;

args_t args;





/*
   This short description is used to generate the output of `bcftools plugin -l`.
   */
const char *about(void)
{
  return
    "Consider only the bi-allelic case. Please use bcftools norm -m-\n";
}

/*
   Called once at startup, allows to initialize local variables.
   Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
   */
int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
  memset(&args,0,sizeof(args_t));
  args.argc = argc; 
  args.argc   = argc; args.argv = argv;
  args.hdr = in;
  printf("pos\tref\talt\ttype\taa\tab\tbb\tnmiss\n");

  return 1;
}

typedef struct
{
    int a, b, phased;
}
gt_t;

int parse_genotype(gt_t *gt, int32_t *ptr);

inline int parse_genotype(gt_t *gt, int32_t *ptr)
{
    if ( ptr[0]==bcf_gt_missing ) return 0;
    if ( ptr[1]==bcf_gt_missing ) return 0;
    if ( ptr[1]==bcf_int32_vector_end ) return 0;
    gt->phased = bcf_gt_is_phased(ptr[1]) ? 1 : 0;
    gt->a = bcf_gt_allele(ptr[0]); if ( gt->a > 1 ) return 0;  // consider only the first two alleles at biallelic sites
    gt->b = bcf_gt_allele(ptr[1]); if ( gt->b > 1 ) return 0;
    return 1;
}


/*
   Called for each VCF record. Return rec to output the line or NULL
   to suppress output.
   */
bcf1_t *process(bcf1_t *rec)
{
    int ngt = bcf_get_genotypes(args.hdr, rec, &args.gt_arr, &args.mgt_arr);
    if ( ngt<0 ) return NULL;
    ngt /= bcf_hdr_nsamples(args.hdr);
    if ( ngt!=2 ) return NULL;
    gt_t gt; 
    int aa, ab, bb, nmiss;
    aa=0; ab=0; bb=0; nmiss=0;
  //output header
  int type = bcf_get_variant_types(rec);
  for (int i =0; i < rec->n_sample; ++i) {
    if(!parse_genotype(&gt, args.gt_arr + i * ngt)){ nmiss++; continue;} 
    if ((gt.a + gt.b) == 0) aa++;
    if ((gt.a + gt.b) == 1) ab++;
    if ((gt.a + gt.b) == 2) bb++;
  }
  printf("%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\n",rec->pos+1, rec->d.allele[0], 
      rec->d.allele[1], type, aa, ab, bb, nmiss); 
  return NULL;
}


/*
   Clean up.
   */
void destroy(void)
{
  free(args.gt_arr);
}


