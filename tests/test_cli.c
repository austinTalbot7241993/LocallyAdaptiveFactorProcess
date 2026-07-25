#include <stdio.h>
#include <stdlib.h>
#include <lafp/gmlib.h>

/**
 * @file test_cli.c
 * @brief Integration test executing the lafp-fit CLI binary.
 */

int main(int argc, char **argv) {
    const char *cli_bin = (argc > 1) ? argv[1] : "../lafp-fit";
    char cmd[1024];

    // Test help command (-h)
    snprintf(cmd, sizeof(cmd), "%s -h > /dev/null", cli_bin);
    if (system(cmd) != 0) GMERR(-10);

    // Test version command (-v)
    snprintf(cmd, sizeof(cmd), "%s -v > /dev/null", cli_bin);
    if (system(cmd) != 0) GMERR(-11);

    // Test fitting with synthetic data files
    snprintf(cmd, sizeof(cmd), "%s -y y.txt -t tobs.txt -n 50 -o test_cli_out > /dev/null", cli_bin);
    if (system(cmd) != 0) GMERR(-20);

    // Verify output file creation
    FILE *fp_theta = fopen("test_cli_out_Theta.txt", "r");
    if (fp_theta == NULL) GMERR(-30);
    fclose(fp_theta);

    FILE *fp_sig = fopen("test_cli_out_Sig.txt", "r");
    if (fp_sig == NULL) GMERR(-31);
    fclose(fp_sig);

    // Remove test outputs
    remove("test_cli_out_Theta.txt");
    remove("test_cli_out_Sig.txt");

    return 0;

GMERRH("main", 1);
}
