# fba; directional
fba map \
    -1 $READ1 \
    -2 $READ2 \
    -w $CELL_BARCODE \
    -f $FEATURE_REF \
    -o matrix_featurecount.csv.gz \
    -r1_coords 0,16 \
    -cb_m 1 \
    -al bwa \
    --mapq 10 \
    -us 16 \
    -ul 12 \
    -um 1 \
    -ud directional \
    --output_directory barcode_mapping


# fba; percentile
fba map \
    -1 $READ1 \
    -2 $READ2 \
    -w $CELL_BARCODE \
    -f $FEATURE_REF \
    -o matrix_featurecount.csv.gz \
    -r1_coords 0,16 \
    -cb_m 1 \
    -al bwa \
    --mapq 10 \
    -us 16 \
    -ul 12 \
    -um 1 \
    -ud percentile \
    --output_directory barcode_mapping


# 10x Cell Ranger
date=`date +%Y%m%d_%H%M%S`
$CELLRANGER count \
    --id=$date \
    --target-panel=Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature_target_panel.csv \
    --transcriptome=$GENOME_REFERENCE_DIR \
    --fastqs=$SEQ_DIR \
    --sample=Targeted_NGSC3_DI_HodgkinsLymphoma_GeneSignature \
    --expect-cells=3394 \
    --nosecondary
