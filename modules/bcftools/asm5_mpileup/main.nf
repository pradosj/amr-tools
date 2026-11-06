nextflow.enable.types = true

// Simple variant calling from a BAM between 2 assembled genomes
process BCFTOOLS_ASM5_MPILEUP {
    container 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'
    memory 10.GB
    cpus 4
    time 1.h

    input:
    record(meta: Record, ref: Path, bam: Path, gff: Path?, args: String?)

    stage:
    stageAs(ref, 'ref.fasta')
    stageAs(bam, 'alignment.bam')
    stageAs(gff, 'ref.gff')

    output:
    record(
        meta: meta,
        vcf: file('mutations.vcf.gz'),
        csi: file('mutations.vcf.gz.csi'),
        txt: file('mutations.txt')
    )

    script:
    def csq_cmd = gff ? "| bcftools csq -Oz --force --local-csq --fasta-ref='ref.fasta' --gff-annot='ref.gff'" : ''
    def csq_qry = gff ? '\\t%INFO/BCSQ' : ''
    """
    bcftools mpileup -Ou \\
        --threads ${task.cpus} \\
        ${args ?: ''} \\
        --per-sample-mF \\
        -a FORMAT/AD \\
        --fasta-ref="ref.fasta" \\
        alignment.bam \\
        | bcftools norm -Ou -m- -f "ref.fasta" \\
        | bcftools filter -Ou -i '(FORMAT/AD[:1] >= 1)' \\
        | bcftools +setGT -Oz -- -t a -n c:1 \\
        ${csq_cmd} \\
        > mutations.vcf.gz
    bcftools index mutations.vcf.gz
    bcftools query -f "%CHROM\\t%POS\\t%REF>%ALT${csq_qry}" mutations.vcf.gz > mutations.txt
    """

    stub:
    """
    touch mutations.vcf.gz mutations.vcf.gz.csi mutations.txt
    """
}
