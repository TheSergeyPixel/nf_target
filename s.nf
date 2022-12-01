#!/usr/bin/env nextflow

nextflow.enable.dsl=2

fastq_ch = Channel.fromFilePairs("${params.rawdir}/${params.samples}", checkIfExists: true, flat:true)
fastq_ch = fastq_ch.map{tuple it[0].split('_')[0], it[1], it[2]}
base_recalibrator_ch = Channel.fromPath("${baseDir}/db/bqsr/*.vcf.gz").collect()


process bwa {
	

	publishDir(path: "${params.outdir}/RawBams", mode: 'copy')
	
	input:
		val(reference)
		tuple val(name), val(r1), val(r2)
		

	output:
		path("*.bam"), emit: bam
		path("*.bai")
		val(name), emit: name

	conda "/root/anaconda3/envs/dna"
	
	script:

	"""
	INDEX="\$((zcat ${r1} | head -n 1 | cut -d ":" -f 1 & \
	zcat ${r1} | head -n 1 | cut -d ":" -f 4) | cat | tr '\\n' '.' | tail -c +2 | head -c -1)"

	bwa mem ${reference} ${r1} ${r2} \
	-t ${task.cpus - 1} -R "@RG\\tID:\${INDEX}\\tSM:${name}\\tPL:Illumina" | \
	samtools sort -o ${name}.bam -
	
	samtools index ${name}.bam
	
	"""
}

process rmdup {

	publishDir(path: "${params.outdir}/RMdup", mode: 'copy', pattern: "*.bam")
	publishDir(path: "${params.outdir}/Metrics", mode: 'copy', pattern: "*.txt")

	input:
		val(gatk)
		file(raw_bam)
		val(name)

	output:
		path("*_rmdup.bam"), emit: bam
		path("*.bai")
		path("*_metrics.txt"), emit: txt
		val(name), emit: name
		val true, emit: ready

	conda '/root/anaconda3/envs/gatk'

	script:
	"""
	python ${gatk} --java-options '-Xmx2G -XX:ParallelGCThreads=1' MarkDuplicates  \
		-I ${raw_bam} \
		-M ${name}_dup_metrics.txt \
		-O ${name}_rmdup.bam \
		--REMOVE_DUPLICATES true \
		--CREATE_INDEX true
	"""
}

process baserecalibrator {

	publishDir(path: "${params.outdir}/Recal_tables", mode: 'copy')

	input:
		val(gatk)
		file(rmdup_bam)
		val(reference)
		tuple val(br1), val(br2), val(br3), val(br4)
		val(bed)
		val(name)
		


	output:
		path("*_recal_data.table"), emit: table
		val(name), emit: name
		val true, emit: ready

	conda '/root/anaconda3/envs/gatk'

	script:
	"""
	python ${gatk} --java-options '-Xmx4G' BaseRecalibrator \
		-I ${rmdup_bam} \
		-R ${reference} \
		--known-sites ${br1} \
		--known-sites ${br2} \
		--known-sites ${br3} \
		--known-sites ${br4} \
		-O ${name}_recal_data.table \
		-L ${bed} -ip 150
	"""

}

process recalibrate {

	publishDir(path: "${params.outdir}/BQSRBAMs", mode: 'copy')

	input:
		val(gatk)
		val(reference)
		val(bed)
		val(name)
		file(bam)
		val(bool)

	output:
		path("*bam"), emit: bam
		path("*.bai")
		val(name), emit: name

	conda '/root/anaconda3/envs/gatk'

	script:

	"""
	python ${gatk} --java-options "-Xms2G -Xmx2G" ApplyBQSR \
		-R ${reference} \
		-I ${bam} \
		--bqsr-recal-file ${params.outdir}/Recal_tables/${name}_recal_data.table \
		-O ${name}BQSR.bam \
		-L ${bed} -ip 150
	"""
}

process haplotypecaller {

	publishDir(path: "${params.outdir}/GVCFs", mode: 'copy')

	input:
		val(bed)
		val(reference)
		val(name)
		val(gatk)
		file(bam) 


	
	output:
		path("*.g.vcf.gz"), emit: gvcf
		file("*.g.vcf.gz.tbi") 
		val(name), emit: name

	conda '/root/anaconda3/envs/gatk'
	script:

	"""
	python ${gatk} --java-options '-Xmx2G' HaplotypeCaller \
		-I ${bam} \
		-R ${reference} \
		-O  ${name}.g.vcf.gz \
		-ERC GVCF --native-pair-hmm-threads 1 \
		-L ${bed} -ip 150
	"""

}


process genotypevcfs {

    publishDir(path: "${params.outdir}/VCFs", mode: 'copy')

    input:
        val(reference)
        val(gatk)
        val(name)
        file(gvcf)

    output:
        path("*.vcf.gz"), emit: vcf
        file("*.tbi")
		val(name), emit: name

	conda '/root/anaconda3/envs/gatk'

    script:
	
    """
    python ${gatk} --java-options '-Xmx2G -XX:ParallelGCThreads=1' GenotypeGVCFs \
        -R ${reference} -V ${baseDir}/results/GVCFs/${gvcf} \
        -O ${name}.vcf.gz \
        -A FragmentLength -A BaseQuality -A ClippingRankSumTest -A MappingQuality
    """
}

process variantfiltration {

	publishDir(path: "${params.outdir}/FilteredVCFs", mode: 'copy')

	input:
		val(reference)
		val(gatk)
		file(vcf)
		val(name)

	output:
		path("*.vcf.gz"), emit: vcf
		file("*.tbi")
		val(name), emit: name

	conda '/root/anaconda3/envs/gatk'

	script:

	"""
	python ${gatk} --java-options '-Xmx2G' VariantFiltration \
		-R ${reference} \
		-V ${baseDir}/results/VCFs/${vcf} \
		-O ${name}F.vcf.gz  \
		--filter-name "FShigher60" --filter-expression "FS > 60.0" \
		--filter-name "SORhigher3" --filter-expression "SOR > 3.0" \
		--filter-name "MQlower40" --filter-expression "MQ < 40.0" \
		--filter-name "MQRankSumlower-12_5" --filter-expression "MQRankSum < -12.5" \
		--filter-name "ReadPosRankSumlower-8" --filter-expression "ReadPosRankSum < -8.0" \
		--filter-name "LowDepth" --filter-expression "DP < 15" \
		--filter-name "LowQUAL"  --filter-expression "QUAL < 23"

	"""

}

process createinterval {

	publishDir(path: "${baseDir}/db/ref", mode: 'copy')

	input:
		val(dict)
		val(gatk)
		val(bed)

	output:
		path("TargetRegions.interval_list"), emit: intervals

	conda '/root/anaconda3/envs/gatk'

	script:
	"""
	python ${gatk} BedToIntervalList \
	-I ${bed} \
	-O TargetRegions.interval_list \
	-SD ${dict}
	"""


}

process collectmetrics {

	publishDir(path: "${params.outdir}/Metrics", mode: 'copy')

	input:
		val(gatk)
		val(reference)
		file(bam)
		file(interval)
		val(name)

	output:
		path("*_insert_size_metrics.txt")
		path("*.pdf")
		path("*_alignment_summary_metrics.txt")
		path("*_HS.txt")
		path("*_per_target.txt")
		val(name), emit: name
		val true, emit: ready



	conda '/root/anaconda3/envs/gatk'

	script:

	"""
	python ${gatk} --java-options '-Xmx2G' CollectInsertSizeMetrics \
		-I ${bam} \
		-O ${name}_insert_size_metrics.txt \
		-R ${reference} -H ${name}.pdf    	
	
	python ${gatk} --java-options '-Xmx2G' CollectAlignmentSummaryMetrics \
		-I ${bam} \
		-O ${name}_alignment_summary_metrics.txt \
		-R ${reference}    	
	
	python ${gatk} --java-options '-Xmx2G' CollectHsMetrics \
		-I ${bam} \
		-R ${reference} \
		-O ${name}_HS.txt \
		-BI  ${interval} -TI ${interval} \
		--NEAR_DISTANCE 75 \
		--PER_TARGET_COVERAGE ${name}_per_target.txt

	"""

}

process vcfstats {

	publishDir(path: "${params.outdir}/Metrics", mode: 'copy')

	input:
		file(vcf)
		val(name)

	output:
		path("*.txt")
		val true, emit: ready

	conda '/root/anaconda3/envs/dna'

	script:
	
	"""
	rtg vcfstats ${vcf} > ${name}_vcf_stats.txt
	"""

	

}

process mergestats {

	publishDir(path: "${params.outdir}/Metrics", mode: 'copy')

	input:
		val(rmdup_ready)
		val(metrics_ready)
		val(vcfstats_ready)

	output:
		path("*.tsv")

	conda '/root/anaconda3/envs/gatk'

	script:

	"""
	python ${baseDir}/bin/nf_merger.py \
		-i ${baseDir}/results/Metrics
	"""


}

workflow {

	bwa(

		params.reference,
		fastq_ch,

	)


	rmdup(

		params.gatk,  
		bwa.out.bam.collect(flat:false).flatten(),
		bwa.out.name.collect(flat:false).flatten()

	)



	baserecalibrator(

		params.gatk, 
		rmdup.out.bam, 
		params.reference, 
		base_recalibrator_ch, 
		params.bed, 
		rmdup.out.name


	)

	recalibrate(

		params.gatk, 
		params.reference, 
		params.bed, 
		rmdup.out.name,  
		rmdup.out.bam,
		baserecalibrator.out.ready

	)

	haplotypecaller(

		params.bed, 
		params.reference, 
		recalibrate.out.name, 
		params.gatk, 
		recalibrate.out.bam
		
	)

	genotypevcfs(

		params.reference,
		params.gatk,
		haplotypecaller.out.name,
		haplotypecaller.out.gvcf

	)

	variantfiltration(

		params.reference,
		params.gatk,
		genotypevcfs.out.vcf,
		genotypevcfs.out.name

	)

	createinterval(

		params.dict,
		params.gatk,
		params.bed

	)

	collectmetrics(

		params.gatk,
		params.reference,
		recalibrate.out.bam,
		createinterval.out,
		recalibrate.out.name

	)

	vcfstats(

		variantfiltration.out.vcf,
		variantfiltration.out.name

	)

	mergestats(

		rmdup.out.ready.collect(),
		collectmetrics.out.ready.collect(),
		vcfstats.out.ready.collect()

	)




	
}

