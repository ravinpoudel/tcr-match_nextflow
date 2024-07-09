# Tcr Match Nextflow

nextflow kuberun https://github.com/ravinpoudel/tcr-match_nextflow.git \
-r main \
-v fsx:/fsx \
--outdir /fsx/pipeline_output/output \
--tcrmatch_database /fsx/pipeline_output/input_data/tcrmatch/refdata/IEDB_data.tsv \
--tcrbetaseq /fsx/pipeline_output/input_data/tcrmatch/tcrseq \
--genomes_dir /fsx/pipeline_output/input_data/tcrmatch/genomes \
-bg -with-trace -with-report -with-dag flowchart.png > log/d1.out & 