#!/usr/bin/env nextflow

import Helper
import CollectInitialMetadata

// Pipeline version
if (workflow.commitId){
    version = "0.1 $workflow.revision"
} else {
    version = "0.1 (local version)"
}

params.help = false
if (params.help){
    Help.print_help(params)
    exit 0
}

def infoMap = [:]
if (params.containsKey("fastq")){
    infoMap.put("fastq", file(params.fastq).size())
}
if (params.containsKey("fasta")){
    if (file(params.fasta) instanceof LinkedList){
        infoMap.put("fasta", file(params.fasta).size())
    } else {
        infoMap.put("fasta", 1)
    }
}

// checks if params.accessions is different from null
if (params.accessions) {
    BufferedReader reader = new BufferedReader(new FileReader(params.accessions));
    int lines = 0;
    while (reader.readLine() != null) lines++;
    reader.close();
    infoMap.put("accessions", lines)

}

Help.start_info(infoMap, "$workflow.start", "$workflow.profile")
CollectInitialMetadata.print_metadata(workflow)


if (!params.accessions && !params.fastq){ exit 1, "'accessions' or 'fastq' parameter missing" }

if (params.accessions){
    IN_accessions_raw = Channel.fromPath(params.accessions).ifEmpty { exit 1, "No accessions file provided with path:'${params.accessions}'" }


    process fasterqDump {

    tag { accession_id }
    publishDir "reads", pattern: "${accession_id}/*fq.gz" // needed?
    maxRetries 1 // adjust?

    input:
    val accession_id from IN_accessions_raw.splitText(){ it.trim() }.filter{ it.trim() != "" }

    output:
    set val({ "$name" != "null" ? "$name" : "$accession_id" }), file("${accession_id}/*fq.gz") optional true into reads_download_out_1_0

    script:
    """
    {
        echo "Downloading the following accession: ${accession_id}"
        fasterq-dump ${accession_id} -e ${task.cpus} -p
        if [ ${params.compress_fastq} = true ]
        then
            echo "Compressing FastQ files..."
            if [ -f ${accession_id}_1.fastq ]
            then
                pigz -p ${task.cpus} ${accession_id}_1.fastq ${accession_id}_2.fastq
            elif [ -f ${accession_id}_3.fastq ]
            then
                echo "No paired end reads were found to compress."
                pigz -p ${task.cpus} ${accession_id}_3.fastq
            else
                echo "FastQ files weren't compressed. Check if FastQ files were downloaded."
            fi
        else
            echo "FastQ files won't be compressed because compress_fastq options was set to: '${params.compress_fastq}.'"
        fi
    } || {
        # If exit code other than 0
        if [ \$? -eq 0 ]
        then
            echo "pass" > .status
        else
            echo "fail" > .status
            echo "Could not download accession $accession_id" > .fail
        fi
    }
    """

    }
}
