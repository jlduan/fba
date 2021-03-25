# fba; percentile
fba extract \
    -1 $READ1 \
    -2 $READ2 \
    -w $CELL_BARCODE \
    -f $FEATURE_BARCODE \
    -o feature_barcoding_output.tsv.gz \
    -r1_coords 0,16 \
    -r2_coords 10,25 \
    -cb_m 2 \
    -fb_m 1

fba count \
    -i feature_barcoding_output.tsv.gz \
    -o matrix_featurecount.csv.gz \
    -us 16 \
    -ul 12 \
    -um 1 \
    -ud percentile


# fba; directional
fba extract \
    -1 $READ1 \
    -2 $READ2 \
    -w $CELL_BARCODE \
    -f $FEATURE_BARCODE \
    -o feature_barcoding_output.tsv.gz \
    -r1_coords 0,16 \
    -r2_coords 10,25 \
    -cb_m 2 \
    -fb_m 1

fba count \
    -i feature_barcoding_output.tsv.gz \
    -o matrix_featurecount.csv.gz \
    -us 16 \
    -ul 12 \
    -um 1 \
    -ud directional


# 10x Cell Ranger
$CELLRANGER count \
    --id=${date} \
    --libraries=samples.csv \
    --feature-ref=feature_references.csv \
    --transcriptome $GENOME_REFERENCE_DIR \
    --expect-cells 1200 \
    --nosecondary


# CITE-seq-Count
CITE-seq-Count \
    -R1 $READ1 \
    -R2 $READ2 \
    -t feature_barcodes_csc.csv \
    -cbf 1 \
    -cbl 16 \
    -umif 17 \
    -umil 28 \
    -cells 1200 \
    -wl cell_barcodes.txt \
    --bc_collapsing_dist 2 \
    --umi_collapsing_dist 1 \
    --max-error 1 \
    -o $date
