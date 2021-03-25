# fba; percentile
fba extract \
    -1 $READ1 \
    -2 $READ2 \
    -w $BARCODES \
    -f $FEATURE_REF \
    -o feature_barcoding_output.tsv.gz \
    -cb_m 1 \
    -fb_m 1 \
    -cb_n 3 \
    -fb_n 3 \
    -r1_coords 0,16 \
    -r2_coords 0,6

fba count \
    -i feature_barcoding_output.tsv.gz \
    -o matrix_featurecount.csv.gz \
    -us 16 \
    -ul 9 \
    -um 1 \
    -ud percentile


# fba; directional
fba extract \
    -1 $READ1 \
    -2 $READ2 \
    -w $BARCODES \
    -f $FEATURE_REF \
    -o feature_barcoding_output.tsv.gz \
    -cb_m 1 \
    -fb_m 1 \
    -cb_n 3 \
    -fb_n 3 \
    -r1_coords 0,16 \
    -r2_coords 0,6

fba count \
    -i feature_barcoding_output.tsv.gz \
    -o matrix_featurecount.csv.gz \
    -us 16 \
    -ul 9 \
    -um 1 \
    -ud directional


# CITE-seq-Count
CITE-seq-Count \
    -R1 $READ1 \
    -R2 $READ2 \
    -t feature_barcodes_csc.csv \
    -cbf 1 \
    -cbl 16 \
    -umif 17 \
    -umil 25 \
    -cells 8617 \
    -wl cell_barcodes.txt \
    --bc_collapsing_dist 1 \
    --umi_collapsing_dist 1 \
    --max-error 1 \
    -o $date
