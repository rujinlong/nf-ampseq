process PREPARE_INPUT {
    publishDir "$baseDir", mode: 'copy', overwrite: true
    publishDir "$params.report", pattern: "metadata.tsv"

    output:
    path("manifest.tsv"), emit: manifest_ch
    path("metadata.tsv"), emit: metadata_ch

    """
    create_metadata.sh $params.raw_seq_dir metadata.csv $params.softlink
    create_metadata.py -i metadata.csv -o metadata.tsv

    echo "sampleid,forward,reverse" > manifest.csv
    for sid in \$(ls $params.raw_seq_dir | sed 's/_R.*//' | sort -u);do
        echo "\$sid,\${sid}_R1${params.raw_seq_ext},\${sid}_R2${params.raw_seq_ext}" >> manifest.csv
    done
    create_manifest.py -i manifest.csv -d $params.raw_seq_dir -o manifest.tsv
    """
}