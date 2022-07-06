#!/bin/bash
#SBATCH --qos=1day
#SBATCH --time=24:00:00
#SBATCH --mem=40g
#SBATCH --output=run.out
#SBATCH --error=run.error
#SBATCH --job-name=demultiplexing
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=jan.waelchli@unibas.ch
#SBATCH --mail-type=ALL

#load modules
module load foss/2018b #interpreters
module load FastQC/0.11.8-Java-1.8
module load cutadapt/2.10-foss-2018b-Python-3.6.6
module load xlsx2csv/0.7.4-foss-2018b-Python-3.6.6

## --------------------------------------------------------------------
## Jan Wälchli | 24.11.2021 | Version 4.0
## --------------------------------------------------------------------

#running time notification
echo 'Start script'

#convert the design file from xlsx to tab
xlsx2csv ../../1_start/design.xlsx design.csv
cat design.csv | tr ',' '\t' | tail -n +2 > design.tab
rm design.csv

#get the runs
 runs=$(awk '{print $3}' design.tab | sort | uniq)

## --------------------------------------------------------------------
## A | Quality Control - FastQC
## --------------------------------------------------------------------

#create output folder
mkdir ../../4_output 2> /dev/null #suppress error message
rm -r  ../../4_output/qc 2> /dev/null
mkdir ../../4_output/qc

#quality control
fastqc -t 20 -k 0 -q ../../2_data/* -o ../../4_output/qc 2> /dev/null

#remove no longer needed files
rm ../../4_output/qc/*.zip

#running time notification
echo 'A - Quality Control done'

## --------------------------------------------------------------------
## B | Barcode Files
## --------------------------------------------------------------------

#create folder
rm -r  demultiplexed 2> /dev/null
mkdir demultiplexed
mkdir demultiplexed/barcodes

#create files
for run in ${runs}; do

	#forward barcodes
	grep ${run} design.tab | awk '{print $6 "-" $8, $7}' | \
	sort | uniq | tr ' ' '\n' | \
	sed 's'/'^F'/'>'${run}'-F'/'g' > demultiplexed/barcodes/${run}_F_barcodes.fasta

	#reverse primers
	grep ${run} design.tab | awk '{print $10 "-" $12, $11}' | \
	sort | uniq | tr ' ' '\n' | \
	sed 's'/'^R'/'>'${run}'-R'/'g' > demultiplexed/barcodes/${run}_R_barcodes.fasta

done

#running time notification
echo 'B - Barcode Files done'

## ---------------------------------------------------------------------
## C | Demultiplexing
## ---------------------------------------------------------------------

for run in ${runs}; do

	mkdir demultiplexed/${run}

	cutadapt \
		-e 0.1 --no-indels \
		-g file:demultiplexed/barcodes/${run}_F_barcodes.fasta \
		-G file:demultiplexed/barcodes/${run}_R_barcodes.fasta \
		-o demultiplexed/${run}/{name1}-{name2}-r1.fastq.gz \
		-p demultiplexed/${run}/{name1}-{name2}-r2.fastq.gz \
		../../2_data/${run}_r1.fastq.gz ../../2_data/${run}_r2.fastq.gz \
		--discard-untrimmed
done

#running time notification
echo 'C - Demultiplexing done'

# ---------------------------------------------------------------------
# D | Clean up
# ---------------------------------------------------------------------

#move samples with primer combinations not in the design file

#forward barcode files
awk '{print $3 "-" $6 "-" $8 "-" $3 "-" $10 "-" $12}' design.tab | \
sed 's/$/-r1.fastq.gz/g' | \
uniq > samples_to_keep.txt

#reverse barcode files
awk '{print $3 "-" $6 "-" $8 "-" $3 "-" $10 "-" $12}' design.tab | \
sed 's/$/-r2.fastq.gz/g' | \
uniq >> samples_to_keep.txt

#move samples with not used barcode combinations
cd demultiplexed
rm -r unsued_samples 2> /dev/null
mkdir unused_samples
for run in ${runs}; do
	cd ${run}
	for sample in *.fastq.gz; do
		if ! grep -qxFe "${sample}" ../../samples_to_keep.txt; then
			mv ${sample} ../unused_samples
		fi
	done
	cd ..
done
cd ..
#not longer required
#rm samples_to_keep.txt

#sort files by taxa

for run in ${runs}; do

	taxa=$(grep ${run} design.tab | awk '{print $2}' | sort | uniq)

	#bacteria
	 if [[ ${taxa} == *'b'* ]]; then
		mkdir demultiplexed/bacteria 2> /dev/null #may already exist
		mkdir demultiplexed/bacteria/${run};
		bac_primers=$(grep ${run} design.tab | awk '{if ($2 == "b") print $0;}' | awk '{print $8}' | sort | uniq)
		for b in ${bac_primers}; do
			mv demultiplexed/${run}/*${b}*.fastq.gz demultiplexed/bacteria/${run}
		done
	fi

	#fungi
	 if [[ ${taxa} == *'f'* ]]; then
		mkdir demultiplexed/fungi 2> /dev/null #may already exist
		mkdir demultiplexed/fungi/${run};
		fun_primers=$(grep ${run} design.tab | awk '{if ($2 == "f") print $0;}' | awk '{print $8}' | sort | uniq)
		for f in ${fun_primers}; do
			mv demultiplexed/${run}/*${f}*.fastq.gz demultiplexed/fungi/${run}
		done
	fi
	rm -r demultiplexed/${run}

done

#running time notification
echo 'D - Clean up done'

# --------------------------------------------------------------------
# E | Primer Files
# --------------------------------------------------------------------

#create folder
rm -r primer_cutted 2> /dev/null
mkdir primer_cutted
mkdir primer_cutted/primers

#create files
for run in ${runs}; do

	taxa=$(grep ${run} design.tab | awk '{print $2}' | sort | uniq)

	#bacteria
	 if [[ ${taxa} == *'b'* ]]; then
		 grep ${run} design.tab | awk '{if ($2 == "b") print $0;}' | \
		 awk '{print $8, $9, $12, $13}' | sort | uniq \
		 > primer_cutted/primers/${run}_b_primers.txt
	 fi

	 #fungi
		if [[ ${taxa} == *'f'* ]]; then
			grep ${run} design.tab | awk '{if ($2 == "f") print $0;}' | \
			awk '{print $8, $9, $12, $13}' | sort | uniq \
			> primer_cutted/primers/${run}_f_primers.txt
		fi

	done

#running time notification
echo 'E - Primer Files done'

# ---------------------------------------------------------------------
# F | Primer cutting
# ---------------------------------------------------------------------

#folders to store cutted sequences
mkdir primer_cutted/bacteria
mkdir primer_cutted/fungi

#repeat for each species in each run
for run in ${runs}; do

	taxa=$(grep ${run} design.tab | awk '{print $2}' | sort | uniq)

		#bacteria
		if [[ ${taxa} == *'b'* ]]; then
			while read p; do #loop over each primer combination
				#get sequence
				f_seq=$(echo ${p} | cut -f2 -d " ")
				r_seq=$(echo ${p} | cut -f4 -d " ")
				#input
				path_in=demultiplexed/bacteria/${run}/
				ls ${path_in} > filenames.txt
				#output
				path_out=primer_cutted/bacteria/${run}/
				mkdir ${path_out}
				#cut primers for each file
				while read infile; do
					#part of the name to keep
					outname=$(echo ${infile} | cut -f1 -d ".")
					direction=$(echo ${infile} | cut -f7 -d "-" | cut -f1 -d ".")
					if [[ ${direction} = "r1" ]]; then
						cutadapt -g ${f_seq} -o ${path_out}${outname}"_cutted".fastq.gz ${path_in}${infile}
					else cutadapt -g ${r_seq} -o ${path_out}${outname}"_cutted".fastq.gz ${path_in}${infile}
					fi
				done < filenames.txt
				rm filenames.txt
			done < primer_cutted/primers/${run}_b_primers.txt
		fi

		#fungi
		if [[ ${taxa} == *'f'* ]]; then
			while read p; do #loop over each primer combination
				#get sequence
				f_seq=$(echo ${p} | cut -f2 -d " ")
				r_seq=$(echo ${p} | cut -f4 -d " ")
				#input
				path_in=demultiplexed/fungi/${run}/
				ls ${path_in} > filenames.txt
				#output
				path_out=primer_cutted/fungi/${run}/
				mkdir ${path_out}
				#cut primers for each file
				while read infile; do
					#part of the name to keep
					outname=$(echo ${infile} | cut -f1 -d ".")
					direction=$(echo ${infile} | cut -f7 -d "-" | cut -f1 -d ".")
					if [[ ${direction} = "r1" ]]; then
						cutadapt -g ${f_seq} -o ${path_out}${outname}"_cutted".fastq.gz ${path_in}${infile}
					else cutadapt -g ${r_seq} -o ${path_out}${outname}"_cutted".fastq.gz ${path_in}${infile}
					fi
				done < filenames.txt
				rm filenames.txt
			done < primer_cutted/primers/${run}_f_primers.txt
		fi

done

echo 'F - primer cutting done'


## ---------------------------------------------------------------------
## F | sequences tracking
## ---------------------------------------------------------------------

rm ../../2_data/seqs.txt 2> /dev/null
touch ../../2_data/seqs.txt

#number of raw sequences per run
for run in ${runs}; do

	#raw
	raw_lines=$(gunzip -c ../../2_data/${run}_r1.fastq.gz | wc -l)
	raw_seqs=$(echo ${raw_lines} / 4 | bc) #each seq has 4 lines

	for species in bacteria fungi; do

		#input
		path_in=primer_cutted/${species}/${run}/

		#count number of demultiplexed sequences
		if [ -d "${path_in}" ]; then #check if path exist, otherwise set to NA
			dem_lines=$(gunzip -c primer_cutted/${species}/${run}/*r1* | wc -l)
			dem_seqs=$(echo ${dem_lines} / 4 | bc) #each seq has 4 lines
		else dem_seqs="NA"
		fi

		#save to corresponding taxa
		if [ ${species} == "bacteria" ]; then
			dem_seqs_bac=$(echo ${dem_seqs})
		else dem_seqs_fun=$(echo ${dem_seqs})
		fi

	done

	#output
	echo ${run} "raw" ${raw_seqs} "bac" ${dem_seqs_bac} "fun" ${dem_seqs_fun} >> ../../2_data/seqs.txt

done

echo 'F - sequences tracking done'
echo 'End script'
