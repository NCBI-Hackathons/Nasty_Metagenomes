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
    set val({ "$name" != "null" ? "$name" : "$accession_id" }), file("${accession_id}/*fq") optional true into IN_fastq_raw

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
        elsenex
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
} else {

   if (!params.fastq){ exit 1, "'fastq' parameter missing"}
   IN_fastq_raw = Channel.fromFilePairs(params.fastq).ifEmpty { exit 1, "No fastq files provided with pattern:'${params.fastq}'" }

}


process fastp {

    tag {sample_id}
    publishDir "results/fastp/", pattern: "*.html"

    input:
    set sample_id, file(fastq_pair) from IN_fastq_raw

    output:
    set sample_id, file('*_QC*')  into OUT_fastq_QC

    script:
    """

    fastp -i ${fastq_pair[0]} -I ${fastq_pair[1]} -o ${sample_id}_QC_1.fq.gz -O ${sample_id}_QC_2.fq.gz -h ${sample_id}.html
    """

}

OUT_fastq_QC.into{ OUT_fastq_QC_1 ; OUT_fastq_QC_2 }

if (params.mode == "magicblast") {

    IN_blastDB = Channel.fromPath("${params.blastdb}*")

    process magicBLAST {

    tag {sample_id}

    input:
    set sample_id, file(fastq_pair) from OUT_fastq_QC_1
    file blastdb from IN_blastDB.collect()

    output:
    set sample_id, file('*_out.sam') into OUT_magicblast

    script:
    """
    echo \$(echo ${blastdb[0]} | sed -r 's/^(.*)\\.[a-z]+\$/\\1/') > blastdb_name.txt
    magicblast -query ${fastq_pair[0]} -query_mate ${fastq_pair[1]} -infmt fastq -db \$(cat blastdb_name.txt) -outfmt sam -out ${sample_id}_out.sam -num_threads ${task.cpus} -paired -no_unaligned
    """
    }

    process samtools {

    tag {sample_id}


    input:
    set sample_id, file(samfile) from OUT_magicblast

    output:
    set sample_id, file('*.depth') into OUT_samtools

    script:
    """
    samtools sort -o sorted.bam -O bam -@ ${task.cpus} ${samfile} && rm *.sam  >> .command.log 2>&1
    samtools depth sorted.bam > ${sample_id}.depth
    """

    }

    IN_reference = Channel.fromPath("${params.reference}")

    process check_coverage {

    tag {sample_id}

    input:
    set sample_id, file(depth_file) from OUT_samtools
    each file(amr_reference) from IN_reference

    output:
    set sample_id, file("*.fasta") into OUT_baits

    script:
    template "magicBlast_depth_parser.pl"
    }

} else if (params.mode == "mash") {

    IN_reference = Channel.fromPath("${params.reference}")

    IN_reference.into{ IN_reference_1 ; IN_reference_2 }

    process mash_sketch {

    tag {amr_reference}
    storeDir 'mash_sketch/'

    input:
    file(amr_reference) from IN_reference_1

    output:
    file("*.msh") into OUT_mash_sketch

    script:
    """
    mash sketch -i ${amr_reference}
    """

    }

    process mash_screen {

    tag {sample_id}

    input:
    set sample_id, file(fastq_pair) from OUT_fastq_QC_1
    file(mash_sketch) from OUT_mash_sketch

    output:
    set sample_id, file("*.screen") into OUT_mash_screen

    script:
    """
    mash screen -p ${task.cpus} -w ${mash_sketch} ${fastq_pair} | awk '\$1>${params.mash_treshold}' > ${sample_id}.amr.screen
    cut -f 5 ${sample_id}.amr.screen
    """
    }

    process mash_hits_compiler {

    tag {sample_id}

    input:
    set sample_id, file(mash_results) from OUT_mash_screen
    file(amr_reference) from IN_reference_2

    output:
    set sample_id, file("*.fasta") into OUT_baits

    script:
    template "mash_hits_compiler.py"
    }



} else if (params.mode == "hmmer") {

    process prepare_fastq {

    tag {sample_id}

    input:
    set sample_id, file(fastq_pair) from OUT_fastq_QC_1

    output:
    set sample_id, file("*.fq") into OUT_preprare_fastq

    script:
    """
    pigz -dc ${fastq_pair} > ${sample_id}.fq
    """
    }

    process read_translator {

    tag {sample_id}

    input:
    set sample_id, file(fastq) from OUT_preprare_fastq

    output:
    set sample_id, file("translated_reads.fa") into OUT_translator

    script:
    template "read_translator.py"
    }

    IN_hmmerDB = Channel.fromPath(params.hmmdb)
    IN_nf_to_seq = Channel.fromPath(params.hmm_fasta)
    IN_threshold = Channel.value(params.hmmr_threshold)

    process hmmer {

    tag {sample_id}

    input:
    set sample_id, file(translated_reads) from OUT_translator
    file hmmerdb from IN_hmmerDB.collect()
    val hmmr_threshold from IN_threshold

    output:
    set sample_id, file("hmmresult*") into OUT_hmmer

    script:
    """
    count=0

    #echo ${hmmerdb}

    for hmmfile in ${hmmerdb};
    do
	    hmmsearch --cpu ${task.cpus} --noali -T ${hmmr_threshold} --acc --incT ${hmmr_threshold} -o "hmmresult_\$count" \$hmmfile ${translated_reads}
	    count=`expr \$count + 1`
    done
    """
    }

    process hits_compiler{

    tag {sample_id}

    input:
    set sample_id, file(hmmresult) from OUT_hmmer
    file nf_to_seq from IN_nf_to_seq

    output:
    set sample_id, file("hit_nucleotide_seqs.fa") into OUT_baits

    script:
    template "hits_compiler.py"
    }


} else {exit 1, "no recognized mode provided. available options: 'magiblast', 'mash', 'hmmer'"}


process guided_assembly {

    tag {sample_id}

    input:
    set sample_id, file(baits), file(fastq_pairs) from OUT_baits.join(OUT_fastq_QC_2)

    output:
    set sample_id, file("*.ga.fa") into OUT_guided_assembly

    script:
    """
    mkdir temp
    cp -L ${fastq_pairs} temp/
    rm ${fastq_pairs} && mv temp/* . && rm -r temp
    gunzip *.gz

    guidedassembler --cores ${task.cpus} --fastq *_1.* --fastq *_2.*  --targets ${baits} --contigs_out ${sample_id}.ga.fa
    """
}

