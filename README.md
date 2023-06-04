# ai-projects
main_assoc
Options:
        -i <INPUT DIR>, --dir=<INPUT DIR>
                path to folder with bed files

        -o <PATH>, --out=<PATH>
                output files directory [default= NULL]

        -g <GENOME.GTF>, --genome=<GENOME.GTF>
                path to genome file

        -c CHROM.SIZES, --chrom=CHROM.SIZES
                path to chrom sizes file

        -v, --verbose
                verbosity level [default= FALSE]

        -p, --simulate
                enable simulated p-value

        -n NUMERIC, --nsim=NUMERIC
                number of simulations [default= 1000]

        -m CHARACTER, --mode=CHARACTER
                set scoring mode/n
              region - use number of overlapping regions /n
              bp - use overlap size /n
              weighted_bp use score column

        -s CHARACTER, --strand=CHARACTER
                ignore strands

        -t CHARACTER, --target=CHARACTER
                path to target file in one vs many mode

        -h, --help
                Show this help message and exit

make_TF_and_hist_table_with_random_ - MC with uniform background

low-input-dif-exp-HSV.R - diffferentian expression analysis
