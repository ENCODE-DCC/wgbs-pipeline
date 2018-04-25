#include "config.h"
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <pthread.h>
#include <signal.h>
#include <locale.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>

#include "libhdr.h"
#include "utils.h"
#include "lk_malloc.h"
#include "loki_compress.h"
#include "bin_tree.h"
#include "string_utils.h"
#include "lkgetopt.h"

#define LOG10 (2.30258509299404568402)
#define DEFAULT_PHRED (0)
#define DEFAULT_REF_PRIOR (1.0)

static void usage(FILE *f) {
  fputs("usage:\n filter_vcf <input file> \n", f);

  fputs("  -o|--out_prefix     PREFIX Prefix name to the output file \n",f);

  fputs("  -a|--all        Output file with all cytosines as well as CpG only "
        "file\n",
        f);
  fprintf(f, "  -t|--threshold  PHRED threshold for genotype calling       "
             "(default='%d')\n",
          DEFAULT_PHRED);
  fprintf(f, "  -r|--ref_prior  Prior weight on 0/0 (reference homozygote) "
             "(default='%g')\n",
          DEFAULT_REF_PRIOR);
  fputs("  -h|help|usage  Print this help \n\n", f);
}

typedef struct {
  char * out_prefix;
  bool all_flag;
  bool recalc_like;
  int threshold;
  double ref_prior;
} filter_params;

filter_params params = {
    .all_flag = false,
    .recalc_like = false,
    .threshold = DEFAULT_PHRED,
    .ref_prior = DEFAULT_REF_PRIOR,
    .out_prefix = "",
};

typedef struct {
  uint64_t x;
  uint64_t counts[8];
  char ref_ctxt[5];
  char alls[3];
  double gl[6];
} vcf_line;

static int cmb_phred[100][100];

static char iupac_cd[256] = {['A'] = 1, ['B'] = 14, ['C'] = 2,  ['D'] = 13,
                             ['G'] = 4, ['H'] = 11, ['K'] = 12, ['M'] = 3,
                             ['R'] = 5, ['S'] = 6,  ['T'] = 8,  ['U'] = 8,
                             ['V'] = 7, ['W'] = 9,  ['Y'] = 10};

static void fill_phred_table(void) {
  for (int i = 0; i < 100; i++) {
    double p1 = exp(-.1 * (double)i * LOG10);
    for (int j = 0; j <= i; j++) {
      double p2 = exp(-.1 * (double)j * LOG10);
      double p = 1.0 - (1.0 - p1) * (1.0 - p2);
      int ph = (int)(0.5 - 10.0 * log(p) / LOG10);
      cmb_phred[i][j] = cmb_phred[j][i] = ph;
    }
  }
}

static int parse_format(char *tp, int *ixp) {
  ixp[0] = ixp[1] = ixp[2] = -1;
  int ix = 0;
  while (*tp) {
    char *tp1 = tp;
    while (*tp1 && *tp1 != ':')
      tp1++;
    if (tp1 - tp == 2) {
      if (tp[0] == 'G' && tp[1] == 'L')
        ixp[0] = ix;
      else if (tp[0] == 'C' && tp[1] == 'X')
        ixp[1] = ix;
    } else if (tp1 - tp == 3 && tp[0] == 'M' && tp[1] == 'C' && tp[2] == '8')
      ixp[2] = ix;
    ix++;
    tp = (*tp1) ? tp1 + 1 : tp1;
  }
  if (ixp[0] < 0 || ixp[1] < 0 || ixp[2] < 0)
    ix = 0;
  return ix;
}

static int get_counts(char *tp, uint64_t *counts) {
  int err = 0;
  char *p;
  int i;
  for (i = 0; i < 7; i++) {
    counts[i] = strtoul(tp, &p, 10);
    if (*p != ',')
      break;
    tp = p + 1;
  }
  if (i == 7) {
    counts[i] = strtoul(tp, &p, 10);
  }
  if (i != 7 || *p)
    err = 1;
  return err;
}

static void output_cpg(FILE *fp, const char *ctg, const uint64_t x,
                       const int phred, const char *ctxt, const char *call_ctxt,
                       const uint64_t *counts, const uint64_t *counts1,
                       const double under_conv, const double over_conv,
                       const int bq_thresh) {
  fprintf(fp, "%s\t%llu\t%.2s\t%.2s\t%d", ctg, x, ctxt + 2, call_ctxt + 2,
          phred);
  uint64_t ct[4];
  if (counts1) {
    ct[0] = counts[5] + counts1[6];
    ct[1] = counts[7] + counts1[4];
    ct[2] = counts[1] + counts[5] + counts[7] + counts1[2] + counts1[4] +
            counts1[7];
    ct[3] = 0;
    for (int i = 0; i < 8; i++)
      ct[3] += counts[i] + counts1[i];
  } else {
    ct[0] = counts[5];
    ct[1] = counts[7];
    ct[2] = counts[1] + counts[5] + counts[7];
    ct[3] = 0;
    for (int i = 0; i < 8; i++)
      ct[3] += counts[i];
  }
  if (ct[0] + ct[1]) {
    double alpha = 1.0 + (double)ct[0];
    double beta = 1.0 + (double)ct[1];
    double m = (alpha - 1.0) / (alpha + beta - 2.0);
    double sd = sqrt(alpha * beta /
                     ((alpha + beta) * (alpha + beta) * (alpha + beta + 1)));
    fprintf(fp, "\t%.3f\t%.3f", m, sd);
  } else
    fputs("\t-\t-", fp);
  fprintf(fp, "\t%llu\t%llu\t%llu\t%llu", ct[0], ct[1], ct[2], ct[3]);
  fprintf(fp, "\t%llu", counts[0]);
  for (int i = 1; i < 8; i++)
    fprintf(fp, ",%llu", counts[i]);
  if (counts1) {
    fprintf(fp, "\t%llu", counts1[0]);
    for (int i = 1; i < 8; i++)
      fprintf(fp, ",%llu", counts1[i]);
  } else {
  	fputs("\t-,-,-,-,-,-,-,-",fp);
  }
  fputc('\n', fp);
}

static void process_file(FILE *fp, filter_params *params) {
  int err = 0;
  tokens *tok = 0, *tok1 = 0, *tok2 = 0, *tok3 = 0;
  string *s = 0, *s1 = 0;
  char *sample = 0;
  void *tbuf = 0;
  FILE *fp_cpg = 0, *fp_all = 0;
  double over_conv = 0.0, under_conv = 0.0;
  int bq_thresh = 0;
  fill_phred_table();
  while (!err) {
    s = fget_string(fp, s, &tbuf);
    if (!s->len)
      break;
    if (!strncmp(get_cstring(s), "##source=bs_call", 16)) {
      tok = tokenize(get_cstring(s), ',', tok);
      for (int i = 1; i < tok->n_tok; i++) {
        if (!strncmp(tok->toks[i], "under_conversion=", 17))
          under_conv = atof(tok->toks[i] + 17);
        else if (!strncmp(tok->toks[i], "over_conversion=", 16))
          over_conv = atof(tok->toks[i] + 16);
        else if (!strncmp(tok->toks[i], "bq_thresh=", 10))
          bq_thresh = atoi(tok->toks[i] + 10);
      }
    } else if (!strncmp(get_cstring(s), "#CHROM", 6)) {
      tok = tokenize(get_cstring(s), '\t', tok);
      if (tok->n_tok >= 10) {
        sample = strdup(tok->toks[9]);
        struct lk_compress *lkc = init_compress();
        char *filter = 0, *suffix = 0;
        if (lkc->default_compress < COMPRESS_NONE) {
          int i = lkc->default_compress;
          filter = lkc->comp_path[i][i == COMPRESS_ZIP ? 1 : 0];
          suffix = lkc->compress_suffix[i];
        }
        if (filter) {
          char *outfile = 0;
          asprintf(&outfile, "%s/%s_cpg.txt.%s",params->out_prefix,sample, suffix);
          int i = child_open(WRITE, outfile, filter);
          fp_cpg = fdopen(i, "w");
          free(outfile);
          if (params->all_flag) {
            asprintf(&outfile, "%s/%s_all.txt.%s",params->out_prefix, sample, suffix);
            i = child_open(WRITE, outfile, filter);
            fp_all = fdopen(i, "w");
            free(outfile);
          }
        } else {
          char *outfile = 0;
          asprintf(&outfile, "%s/%s_cpg.txt",params->out_prefix,sample);
          fp_cpg = fopen(outfile, "w");
          free(outfile);
          if (params->all_flag) {
            asprintf(&outfile, "%s/%s_all.txt",params->out_prefix,sample);
            fp_all = fopen(outfile, "w");
            free(outfile);
          }
        }
        break;
      }
    }
  }
  if (!s || !s->len)
    return;
  bool prev_flag = false;
  while (!err) {
    if (!prev_flag) {
      s = fget_string(fp, s, &tbuf);
      if (!s->len)
        break;
      tok = tokenize(get_cstring(s), '\t', tok);
      if (tok->n_tok < 10)
        continue;
    } else
      prev_flag = false;
    char *ctxt = tok->toks[7];
    bool ref_cg = false;
    if (strlen(ctxt) == 8 && !strncmp(ctxt, "CX=", 3)) {
      ctxt += 3;
      if (ctxt[2] == 'C' && ctxt[3] == 'G')
        ref_cg = true;
      int ixp[3];
      int ix = parse_format(tok->toks[8], ixp);
      if (!ix) {
        fprintf(stderr, "Bad format field\n");
        break;
      }
      tok1 = tokenize(tok->toks[9], ':', tok1);
      if (tok1->n_tok != ix) {
        fprintf(stderr, "Bad number of columns in genotype field\n");
        break;
      }
      char *call_ctxt = tok1->toks[ixp[1]];
      if (strlen(call_ctxt) != 5) {
        fprintf(stderr, "Bad CX sub field in genotype field\n");
        break;
      }
      bool call_cg = false;
      if ((iupac_cd[(int)call_ctxt[2]] & 2) &&
          (iupac_cd[(int)call_ctxt[3]] & 4))
        call_cg = true;
      uint64_t counts[8];
      if (get_counts(tok1->toks[ixp[2]], counts)) {
        fprintf(stderr, "Bad format for counts field\n");
        break;
      }
      if (ref_cg || call_cg) {
        int phred = atoi(tok->toks[5]);
        if (phred > 99)
          phred = 99;
        s1 = fget_string(fp, s1, &tbuf);
        if (!s1->len)
          break;
        tok2 = tokenize(get_cstring(s1), '\t', tok2);
        if (tok2->n_tok < 10)
          continue;
        uint64_t x = atol(tok->toks[1]);
        uint64_t x1 = atol(tok2->toks[1]);
        if (!strcmp(tok->toks[0], tok2->toks[0]) && x1 == x + 1) {
          int ixp1[3];
          int ix1 = parse_format(tok2->toks[8], ixp1);
          if (!ix1) {
            fprintf(stderr, "Bad format field\n");
            break;
          }
          tok3 = tokenize(tok2->toks[9], ':', tok3);
          if (tok3->n_tok != ix1) {
            fprintf(stderr, "Bad number of columns in genotype field\n");
            break;
          }
          uint64_t counts1[8];
          if (get_counts(tok3->toks[ixp1[2]], counts1)) {
            fprintf(stderr, "Bad format for counts field\n");
            break;
          }
          int phred1 = atoi(tok2->toks[5]);
          if (phred1 > 99)
            phred1 = 99;
          output_cpg(fp_cpg, tok->toks[0], x, cmb_phred[phred][phred1], ctxt,
                     call_ctxt, counts, counts1, under_conv, over_conv,
                     bq_thresh);
        } else {
          output_cpg(fp_cpg, tok->toks[0], x, phred, ctxt, call_ctxt, counts, 0,
                     under_conv, over_conv, bq_thresh);
          prev_flag = true;
          string *ts = s;
          s = s1;
          s1 = ts;
          tokens *ttok = tok;
          tok = tok2;
          tok2 = ttok;
        }
      }
    }
  }
}

int main(int argc, char *argv[]) {
  static struct option longopts[] = {{"out_prefix", required_argument, 0, 'o'},
                                     {"threshold", required_argument, 0, 't'},
                                     {"ref_prior", required_argument, 0, 'r'},
                                     {"all", no_argument, 0, 'a'},
                                     {"help", no_argument, 0, 'h'},
                                     {"usage", no_argument, 0, 'h'},
                                     {0, 0, 0, 0}};
  int err = 0;
  int c;
  while (!err && (c = getopt_long(argc, argv, "o:t:r:ah?", longopts, 0)) != -1) {
    switch (c) {
    case 'o':
      params.out_prefix = optarg;
      break;
    case 't':
      params.threshold = atoi(optarg);
      break;
    case 'r':
      params.ref_prior = atof(optarg);
      params.recalc_like = true;
      break;
    case 'a':
      params.all_flag = true;
      break;
    case 'h':
    case '?':
      err = 1;
      usage(stdout);
      break;
    }
  }
  if (err == 1)
    return 0;
  FILE *fp;
  if (argc > optind) {
    int j;
    fp = open_readfile_and_check(argv[optind], &j, init_compress());
  } else
    fp = stdin;
  if (fp)
    process_file(fp, &params);
  if (fp != stdin)
    fclose(fp);
  return err;
}
