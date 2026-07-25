#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <lafp/gmlib.h>
#include <lafp/matrix_utils.h>
#include <lafp/mdarray.h>
#include <lafp/NGPmcmc.h>

#define LAFP_VERSION "1.0.0"

static void print_usage(const char *prog_name) {
    printf("Usage: %s -y <y_file> -t <tobs_file> [options]\n\n", prog_name);
    printf("Locally Adaptive Factor Process (LAFP) Non-Gaussian Process MCMC Sampler\n\n");
    printf("Required Arguments:\n");
    printf("  -y, --input-y <path>     Path to observation matrix input file (y.txt)\n");
    printf("  -t, --input-t <path>     Path to observation timestamps input file (tobs.txt)\n\n");
    printf("Sampling & Output Options:\n");
    printf("  -n, --niter <int>        Total MCMC sampling iterations (default: 2000)\n");
    printf("  -o, --out-prefix <str>   Output filename prefix (default: 'lafp_out')\n");
    printf("      --nt <int>           Row count of time points (auto-detected if omitted)\n\n");
    printf("Hyperparameter & Prior Options:\n");
    printf("      --sig-u <val>        Initial state process noise scale sigU (default: 1000.0)\n");
    printf("      --sig-a <val>        Initial derivative process noise scale sigA (default: 5.0)\n");
    printf("      --sig-eps <val>      Initial observation noise std dev sigEps (default: 2.0)\n");
    printf("      --sig-mu <val>       Prior std dev for mean state components sigMu (default: 4.0)\n");
    printf("      --sig-alpha <val>    Prior std dev for derivative state components sigAlph (default: 4.0)\n");
    printf("  -a, --prior-a <val>      Inverse-Gamma prior shape parameter a (default: 1.0)\n");
    printf("  -b, --prior-b <val>      Inverse-Gamma prior scale parameter b (default: 1.0)\n\n");
    printf("General Options:\n");
    printf("  -h, --help               Display this help message and exit\n");
    printf("  -v, --version            Display version information and exit\n");
}

static int count_file_lines(const char *filepath) {
    FILE *fp = fopen(filepath, "r");
    if (fp == NULL) return -1;
    int count = 0;
    int ch;
    while ((ch = fgetc(fp)) != EOF) {
        if (ch == '\n') count++;
    }
    fclose(fp);
    return count;
}

int main(int argc, char **argv) {
    char *y_path = NULL;
    char *tobs_path = NULL;
    char *out_prefix = "lafp_out";
    int Niter = 2000;
    int Nt = -1;

    double sigU_val = 1000.0;
    double sigA_val = 5.0;
    double sigEps_val = 2.0;
    double sigMu_val = 4.0;
    double sigAlph_val = 4.0;
    double a_val = 1.0;
    double b_val = 1.0;

    static struct option long_options[] = {
        {"input-y",     required_argument, 0, 'y'},
        {"input-t",     required_argument, 0, 't'},
        {"niter",       required_argument, 0, 'n'},
        {"out-prefix",  required_argument, 0, 'o'},
        {"nt",          required_argument, 0, 1001},
        {"sig-u",       required_argument, 0, 1002},
        {"sig-a",       required_argument, 0, 1003},
        {"sig-eps",     required_argument, 0, 1004},
        {"sig-mu",      required_argument, 0, 1005},
        {"sig-alpha",   required_argument, 0, 1006},
        {"prior-a",     required_argument, 0, 'a'},
        {"prior-b",     required_argument, 0, 'b'},
        {"help",        no_argument,       0, 'h'},
        {"version",     no_argument,       0, 'v'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, "y:t:n:o:a:b:hv", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'y':
                y_path = optarg;
                break;
            case 't':
                tobs_path = optarg;
                break;
            case 'n':
                Niter = atoi(optarg);
                break;
            case 'o':
                out_prefix = optarg;
                break;
            case 'a':
                a_val = atof(optarg);
                break;
            case 'b':
                b_val = atof(optarg);
                break;
            case 1001:
                Nt = atoi(optarg);
                break;
            case 1002:
                sigU_val = atof(optarg);
                break;
            case 1003:
                sigA_val = atof(optarg);
                break;
            case 1004:
                sigEps_val = atof(optarg);
                break;
            case 1005:
                sigMu_val = atof(optarg);
                break;
            case 1006:
                sigAlph_val = atof(optarg);
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            case 'v':
                printf("lafp-fit version %s\n", LAFP_VERSION);
                return 0;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    if (y_path == NULL || tobs_path == NULL) {
        fprintf(stderr, "Error: Both -y/--input-y and -t/--input-t arguments are required.\n\n");
        print_usage(argv[0]);
        return 1;
    }

    // Auto-detect Nt if omitted
    if (Nt <= 0) {
        Nt = count_file_lines(tobs_path);
        if (Nt <= 0) {
            fprintf(stderr, "Error: Could not read time series line count from '%s'. Please specify --nt.\n", tobs_path);
            return 1;
        }
    }

    printf("[LAFP-FIT] Starting MCMC Fit...\n");
    printf("  Input y      : %s\n", y_path);
    printf("  Input tobs   : %s\n", tobs_path);
    printf("  Time Points  : %d\n", Nt);
    printf("  Iterations   : %d\n", Niter);
    printf("  Out Prefix   : %s\n", out_prefix);

    FILE *fp_y = fopen(y_path, "r");
    if (fp_y == NULL) {
        fprintf(stderr, "Error: Unable to open input file '%s'\n", y_path);
        return 1;
    }
    gsl_matrix *y = gsl_matrix_alloc(Nt, 1);
    if (gsl_matrix_fscanf(fp_y, y)) {
        fprintf(stderr, "Error: Failed to read matrix data from '%s'\n", y_path);
        fclose(fp_y);
        gsl_matrix_free(y);
        return 1;
    }
    fclose(fp_y);

    FILE *fp_t = fopen(tobs_path, "r");
    if (fp_t == NULL) {
        fprintf(stderr, "Error: Unable to open input file '%s'\n", tobs_path);
        gsl_matrix_free(y);
        return 1;
    }
    gsl_vector *tobs = gsl_vector_alloc(Nt);
    if (gsl_vector_fscanf(fp_t, tobs)) {
        fprintf(stderr, "Error: Failed to read vector data from '%s'\n", tobs_path);
        fclose(fp_t);
        gsl_matrix_free(y);
        gsl_vector_free(tobs);
        return 1;
    }
    fclose(fp_t);

    gsl_vector *sigU = gsl_vector_alloc(1);
    gsl_vector *sigA = gsl_vector_alloc(1);
    gsl_vector_set(sigU, 0, sigU_val);
    gsl_vector_set(sigA, 0, sigA_val);

    marray3d *th = marray3d_alloc(Niter, Nt, 3);
    gsl_matrix *sig = gsl_matrix_alloc(Niter, 3);

    NGPmcmc *mcmc = NGPmcmc_New();
    if (mcmc == NULL) {
        fprintf(stderr, "Error: Failed to allocate NGPmcmc sampler\n");
        return 1;
    }

    if (NGPmcmc_construct(mcmc, y, tobs, Niter, sigU, sigA,
                          sigEps_val, sigMu_val, sigAlph_val,
                          a_val, b_val, th, sig)) {
        fprintf(stderr, "Error: Failed to construct NGPmcmc sampler\n");
        NGPmcmc_free(mcmc);
        return 1;
    }

    printf("[LAFP-FIT] Executing MCMC sampling iterations...\n");
    if (NGPmcmc_operations(mcmc)) {
        fprintf(stderr, "Error: MCMC sampling failed during execution\n");
        NGPmcmc_free(mcmc);
        return 1;
    }

    printf("[LAFP-FIT] Exporting output files with prefix '%s'...\n", out_prefix);
    if (NGPmcmc_writeOutputs(mcmc, out_prefix)) {
        fprintf(stderr, "Error: Failed to write output files\n");
        NGPmcmc_free(mcmc);
        return 1;
    }

    printf("[LAFP-FIT] MCMC Fit Completed Successfully.\n");

    NGPmcmc_free(mcmc);
    gsl_matrix_free(y);
    gsl_vector_free(tobs);
    gsl_vector_free(sigU);
    gsl_vector_free(sigA);
    marray3d_free(th);
    gsl_matrix_free(sig);

    return 0;
}
