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
#include <lafp/NGPlaf2.h>

#define LAFP_VERSION "1.0.0"

static void print_usage(const char *prog_name) {
    printf("Usage: %s -y <y_file> -t <tobs_file> [options]\n\n", prog_name);
    printf("Locally Adaptive Factor Process (LAFP) Sampler Engine\n");
    printf("Supports both Multivariate LAFP (Durante et al., JMLR 2014) and Univariate nGP (Zhu & Dunson, JASA 2013).\n\n");
    printf("Required Arguments:\n");
    printf("  -y, --input-y <path>     Path to observation matrix input file (y.txt)\n");
    printf("  -t, --input-t <path>     Path to observation timestamps input file (tobs.txt)\n\n");
    printf("Sampling & Output Options:\n");
    printf("  -n, --niter <int>        Total MCMC sampling iterations (default: 2000)\n");
    printf("  -o, --out-prefix <str>   Output filename prefix (default: 'lafp_out')\n");
    printf("      --nt <int>           Row count of time points (auto-detected if omitted)\n");
    printf("  -p, --np <int>           Column count of response variables (auto-detected if omitted)\n\n");
    printf("Multivariate Model Options (Durante et al. 2014):\n");
    printf("  -k, --n-factors <int>    Number of latent factors K (default: auto for p > 1)\n");
    printf("  -l, --n-dict <int>       Number of dictionary process elements L (default: 2)\n\n");
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

static int count_file_cols(const char *filepath) {
    FILE *fp = fopen(filepath, "r");
    if (fp == NULL) return -1;
    char line[65536];
    int cols = 0;
    while (fgets(line, sizeof(line), fp)) {
        char *ptr = line;
        while (*ptr == ' ' || *ptr == '\t' || *ptr == '\r' || *ptr == '\n') ptr++;
        if (*ptr == '\0') continue;
        int in_token = 0;
        while (*ptr && *ptr != '\r' && *ptr != '\n') {
            if (*ptr == ' ' || *ptr == '\t') {
                in_token = 0;
            } else if (!in_token) {
                in_token = 1;
                cols++;
            }
            ptr++;
        }
        break;
    }
    fclose(fp);
    return cols;
}

int main(int argc, char **argv) {
    char *y_path = NULL;
    char *tobs_path = NULL;
    char *out_prefix = "lafp_out";
    int Niter = 2000;
    int Nt = -1;
    int Np = -1;
    int NK = 0;
    int NL = 2;

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
        {"np",          required_argument, 0, 'p'},
        {"n-factors",   required_argument, 0, 'k'},
        {"n-dict",      required_argument, 0, 'l'},
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

    while ((opt = getopt_long(argc, argv, "y:t:n:o:p:k:l:a:b:hv", long_options, &option_index)) != -1) {
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
            case 'p':
                Np = atoi(optarg);
                break;
            case 'k':
                NK = atoi(optarg);
                break;
            case 'l':
                NL = atoi(optarg);
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

    // Auto-detect Nt and Np
    if (Nt <= 0) {
        Nt = count_file_lines(tobs_path);
        if (Nt <= 0) {
            fprintf(stderr, "Error: Could not read time series line count from '%s'. Please specify --nt.\n", tobs_path);
            return 1;
        }
    }
    if (Np <= 0) {
        Np = count_file_cols(y_path);
        if (Np <= 0) {
            fprintf(stderr, "Error: Could not detect column count from '%s'. Please specify --np.\n", y_path);
            return 1;
        }
    }

    // Decide multivariate (Durante 2014) vs univariate (Zhu 2013)
    int use_multivariate = (Np > 1 || NK > 0);
    if (use_multivariate && NK <= 0) {
        NK = (Np > 1) ? (Np < 2 ? Np : 2) : 1;
    }

    printf("[LAFP-FIT] Starting MCMC Fit...\n");
    printf("  Input y      : %s (%d variables x %d time points)\n", y_path, Np, Nt);
    printf("  Input tobs   : %s (%d time points)\n", tobs_path, Nt);
    printf("  Model Engine : %s\n", use_multivariate ? "Multivariate LAFP (Durante et al. 2014)" : "Univariate nGP (Zhu & Dunson 2013)");
    if (use_multivariate) {
        printf("  Factors (K)  : %d\n", NK);
        printf("  Dict Dim (L) : %d\n", NL);
    }
    printf("  Iterations   : %d\n", Niter);
    printf("  Out Prefix   : %s\n", out_prefix);

    FILE *fp_y = fopen(y_path, "r");
    if (fp_y == NULL) {
        fprintf(stderr, "Error: Unable to open input file '%s'\n", y_path);
        return 1;
    }
    gsl_matrix *y = gsl_matrix_alloc(Nt, Np);
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

    if (use_multivariate) {
        int Nsamps = Niter / 5;
        if (Nsamps < 1) {
            fprintf(stderr, "Error: Niter must be at least 5 for NGPlaf2 sampler.\n");
            gsl_matrix_free(y);
            gsl_vector_free(tobs);
            return 1;
        }

        marray3d *Ksi_out = marray3d_alloc(Nsamps, Nt, NK * NL);
        marray3d *Psi_out = marray3d_alloc(Nsamps, Nt, NK);
        marray3d *yhat_out = marray3d_alloc(Nsamps, Nt, Np);
        marray3d *mu_out = marray3d_alloc(Nsamps, Nt, Np);
        marray4d *Sigma_out = marray4d_alloc(Nsamps, Nt, Np, Np);

        gsl_vector *sigPrior = gsl_vector_alloc(2);
        gsl_vector *epsPrior = gsl_vector_alloc(2);
        gsl_vector *ksiPrior = gsl_vector_alloc(2);
        gsl_vector *APrior   = gsl_vector_alloc(2);
        gsl_vector *psiPrior = gsl_vector_alloc(2);
        gsl_vector *BPrior   = gsl_vector_alloc(2);
        gsl_vector *aaPrior  = gsl_vector_alloc(2);

        gsl_vector_set(sigPrior, 0, sigMu_val);
        gsl_vector_set(sigPrior, 1, sigAlph_val);
        gsl_vector_set(epsPrior, 0, 1.0);
        gsl_vector_set(epsPrior, 1, 0.1);
        gsl_vector_set(ksiPrior, 0, a_val);
        gsl_vector_set(ksiPrior, 1, b_val);
        gsl_vector_set(APrior,   0, 1.0);
        gsl_vector_set(APrior,   1, 1.0);
        gsl_vector_set(psiPrior, 0, 1.0);
        gsl_vector_set(psiPrior, 1, 1.0);
        gsl_vector_set(BPrior,   0, 1.0);
        gsl_vector_set(BPrior,   1, 1.0);
        gsl_vector_set(aaPrior,  0, 3.0);
        gsl_vector_set(aaPrior,  1, 3.0);

        NGPlaf2 *laf = NGPlaf2_New();
        if (laf == NULL) {
            fprintf(stderr, "Error: Failed to allocate NGPlaf2 sampler\n");
            return 1;
        }

        if (NGPlaf2_construct(laf, tobs, y, NK, NL, Niter,
                              Ksi_out, Psi_out, yhat_out, mu_out, Sigma_out,
                              sigPrior, epsPrior, ksiPrior, APrior,
                              psiPrior, BPrior, aaPrior)) {
            fprintf(stderr, "Error: Failed to construct NGPlaf2 sampler\n");
            NGPlaf2_free(laf);
            return 1;
        }

        printf("[LAFP-FIT] Executing Multivariate MCMC sampling iterations...\n");
        if (NGPlaf2_operations(laf)) {
            fprintf(stderr, "Error: Multivariate MCMC sampling failed during execution\n");
            NGPlaf2_free(laf);
            return 1;
        }

        printf("[LAFP-FIT] Exporting output files with prefix '%s'...\n", out_prefix);
        if (NGPlaf2_writeOutputs(laf, out_prefix)) {
            fprintf(stderr, "Error: Failed to write output files\n");
            NGPlaf2_free(laf);
            return 1;
        }

        NGPlaf2_free(laf);
        gsl_vector_free(sigPrior);
        gsl_vector_free(epsPrior);
        gsl_vector_free(ksiPrior);
        gsl_vector_free(APrior);
        gsl_vector_free(psiPrior);
        gsl_vector_free(BPrior);
        gsl_vector_free(aaPrior);
        marray3d_free(Ksi_out);
        marray3d_free(Psi_out);
        marray3d_free(yhat_out);
        marray3d_free(mu_out);
        marray4d_free(Sigma_out);
    } else {
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

        printf("[LAFP-FIT] Executing Univariate MCMC sampling iterations...\n");
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

        NGPmcmc_free(mcmc);
        gsl_vector_free(sigU);
        gsl_vector_free(sigA);
        marray3d_free(th);
        gsl_matrix_free(sig);
    }

    gsl_matrix_free(y);
    gsl_vector_free(tobs);

    printf("[LAFP-FIT] MCMC Fit Completed Successfully.\n");
    return 0;
}
