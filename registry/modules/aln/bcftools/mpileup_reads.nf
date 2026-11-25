nextflow.enable.types = true

// Variant calling (haploid) from a BAM of reads (ont-sup preset by default) aligned on a reference
process BCFTOOLS_MPILEUP_READS {
    container 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'
    memory 10.GB
    cpus 4
    time 2.h

    input:
    record(
        meta: Record,
        ref: Path,          // the reference FASTA the reads are aligned on
        bam: Path,          // the sorted BAM of the reads
        gff: Path?,         // annotation of the reference: consequences of the variants with bcftools csq
        filter: String?,    // bcftools filter expression of the variants kept
        preset: String?,    // passed to -X flag of mpileup, (default: ont-sup)
        args: String?       // other options of bcftools mpileup
    )

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
    def default_args = '--ff=UNMAP --ambig-reads=drop -Q5 --max-BQ 30 --gap-frac=0.01 --min-ireads=5' //--no-BAQ  --indel-size 1000 
    def keep = filter ?: '(FORMAT/ADF[:1] >= DP*0.1) && (FORMAT/ADR[:1] >= DP*0.1)'
    def csq_cmd = gff ? "| bcftools csq -Oz --force --local-csq --fasta-ref='ref.fasta' --gff-annot='ref.gff'" : ''
    def csq_qry = gff ? '\\t%INFO/BCSQ' : ''
    
    """
    bcftools mpileup -Ou \\
        --threads ${task.cpus} \\
        ${args ?: default_args} \\
        -X ${preset ?: 'ont-sup'}
        -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP \\
        --fasta-ref="ref.fasta" \\
        --per-sample-mF \\
        alignment.bam \\
        | bcftools call -Ou -m -v --ploidy 1 \\
        | bcftools norm -Ou -m- -f "ref.fasta" \\
        | bcftools filter -Oz -i '${keep}' \\
        ${csq_cmd} \\
        > mutations.vcf.gz
    bcftools index mutations.vcf.gz
    bcftools query -f "%CHROM\\t%POS\\t%REF>%ALT\\t%QUAL\\t[%DP\\t%AD]${csq_qry}\\n" mutations.vcf.gz > mutations.txt
    """

    stub:
    """
    touch mutations.vcf.gz mutations.vcf.gz.csi mutations.txt
    """
}






