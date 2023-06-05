


# Note: make sure reference files ref_file.gff, ref_file.fasta, snpEff.jar and snpEff.config are in the working directory
# SnpEff version SnpEff 4.3s (build 2017-10-25 10:05), by Pablo Cingolani

ref=$1;
ref_file=$2;
pair_ids=$3;
vcf_dir=$4;

wd="/project/directory/snippy/";
cd $wd

if [ ! -f $wd"/data/"$ref_file"/snpEffectPredictor.bin" ]
then
	# mkdir $wd"data"
	mkdir $wd"data/"$ref_file
	cp $wd$ref_file".gff" $wd"data/"$ref_file"/genes.gff"
	cp $wd$ref_file".fasta" $wd"data/"$ref_file"/sequences.fa"
	java -jar snpEff.jar build -gff3 -v $ref_file
fi

config_filename=$wd'snpEff.config'; # A new line: 'ref_file.genome : ref_file.fa' needs to be added to snpEff.config, where ref_file is the reference file name
genome=$ref_file;

laPairs=`awk -F'\t' '{ print $1}' $pair_ids`;

for pair in $laPairs
do
	echo $pair
	vcf_file=$vcf_dir$pair".snippy."$ref".vcf";
	annotation_stats_file=$vcf_dir$pair".snippy."$ref".snpEff_stats.html";
	ann_vcf=$vcf_dir$pair".snippy."$ref".snpEff_ann.vcf";

	if [ ! -f $ann_vcf ]
	then
		bsub -q normal -G team81 -J $pair"_ann" -o $pair"_ann.out" -R "select[mem > 5000] rusage[mem=5000]" -M 5000 "java -jar $wd'snpEff.jar' ann -nodownload -no-downstream -no-upstream -verbose -stats $annotation_stats_file -c $config_filename $genome $vcf_file > $ann_vcf"
	fi
done