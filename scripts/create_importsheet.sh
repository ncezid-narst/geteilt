#!/bin/bash

#Post-geteilt bash script to compile import sheet for NARMS database update

#HELP function
function HELP {
echo ""
echo "Usage:" $0
echo "			-i /path/to/geteilt/output/"
echo ""
echo "Example: $0 -i /path/to/geteilt/output/"
echo ""
exit 0
}

###Take arguments
#Run HELP if -h -? or invalid input
#Set ASSEMBLY to -a
#Set LONGREADS to -r
while getopts ":hi:" option; do
	case ${option} in
		h)
		HELP
		;;
		i)
		export OUTPUT=${OPTARG}
		;;
		\?)
		echo "Invalid option: ${OPTARG}" 1>&2
		HELP
		;;
	esac
done

#create time stamp
time_stamp=$(date +%Y_%m_%d_%H_%M_%S)

#Get fullpath to geteilt output
FULLOUTPUT=$(realpath ${OUTPUT})

#For each subdirectory under geteilt output directory, do these things:
for subdir in ${FULLOUTPUT}/*/
do
	sample_name=$(basename ${subdir})
	#echo $subdir
	#Create import.resfinder.tsv from resfinder.tsv - first two columns; replace 'contigs' with name of subdirectory; strip white spaces
	#Create staramr.pointfinder.tsv if it's present (only present in Salmonella and Campylobacter screens) - same deal as above
	#Create import.plasmidfinder.tsv from plasmidfinder.tsv - same deal as above
	#NOTE 2023/06/12 - PUTTING IN TEMPORARY FIX FOR SORTING DUPLICATE HITS OUT OF RESFINDER.TSV: CHECKS 7TH COLUMN (CONTIG NAME) AND FILTERS BY UNIQUE VALUES
	if [ -d ${subdir}/staramr/ ]; then
		cat ${subdir}/staramr/resfinder.tsv | cut -f 1,2,7 | sed -e "s/contigs/$sample_name/" | sed -e 's/\ //' | sed '1d' | sort -u -k3,3 | cut -f 1-2 > ${subdir}/import.resfinder.tsv
		#Pull "None" from summary.tsv in case resfinder.tsv is empty (can't import null)
		cat ${subdir}/staramr/summary.tsv | awk -F "\t" 'BEGIN {OFS="\t"}{if ($3 == "None") print $1,$3 }' | sed -e "s/contigs/$sample_name/" >> ${subdir}/import.resfinder.tsv
		if [ -e ${subdir}/staramr/pointfinder.tsv ]; then
			cat ${subdir}/staramr/pointfinder.tsv | cut -f 1-2 | sed -e "s/contigs/$sample_name/" | sed -e 's/\ //' | sed '1d'| sed -e 's/([A-Z]/(/' |sed -e 's/[A-Z])/)/' > ${subdir}/staramr.pointfinder.tsv
		else
			echo "boop"
		fi
		cat ${subdir}/staramr/plasmidfinder.tsv | cut -f 1-2 |  sed -e "s/contigs/$sample_name/" | sed -e 's/\ //' | sed '1d' > ${subdir}/import.plasmidfinder.tsv
		#Pull "None" from summary.tsv in case plasmidfinder.tsv is empty (can't import null)
		cat ${subdir}/staramr/summary.tsv | awk -F "\t" 'BEGIN {OFS="\t"}{if ($5 == "None") print $1,$5 }' | sed -e "s/contigs/$sample_name/" >> ${subdir}/import.plasmidfinder.tsv
		cut -f 1,3  ${subdir}/staramr/summary.tsv | sed '1d' | sed -e "s/contigs/$sample_name/" | sed -e 's/\ //g' > ${subdir}/staramr_summary.tsv
	else
		echo "moooooooooooooo"
	fi
	#Create import.mutational.tsv for non-Salmonella
	#Create import.23detection.tsv for Campylobacter coli, Campylobacter jejuni, Campylobacter upsaliensis
	#Create import.gyrAdetection.tsv for Campylobacter coli, Campylobacter upsaliensis
	if [ -d ${subdir}/mutational/ ]; then
		if [ -e ${subdir}/mutational/report.tsv ]; then
			cat ${subdir}/mutational/report.tsv | awk -F "\t" 'BEGIN{OFS="\t"} {if ($18 == "1") {print "sample_name",$1"("$19")"}}' | sed -e "s/sample_name/$sample_name/" | sed -e 's/(./(/' | sed -e 's/.)/)/' > ${subdir}/import.mutational.tsv
			cat ${subdir}/mutational/report.tsv | awk -F "\t" 'BEGIN{OFS="\t"} {if ($18 == "1") {print "sample_name",$1"("$19")"}}' | sed -e "s/sample_name/$sample_name/" > ${subdir}/${sample_name}.summary.mutational.tsv
			if [ -s ${subdir}/${sample_name}.summary.mutational.tsv ]; then
					cat ${subdir}/${sample_name}.summary.mutational.tsv | echo -e "${sample_name}\t""$(cut -f 2 | paste -d, -s)" > ${subdir}/all.summary.mutational.tsv
				else
					echo "heE-hAw"
			fi
			if [ ! -s ${subdir}/${sample_name}.summary.mutational.tsv ]; then
				if [ -f ${subdir}/mutational/log.clusters.gz ]; then
						echo -e "${sample_name}\tNone" > ${subdir}/all.summary.mutational.tsv
						echo -e "${sample_name}\tNone" > ${subdir}/import.mutational.tsv
					else
						echo -e "${sample_name}\tError" > ${subdir}/all.summary.mutational.tsv  
						echo -e "${sample_name}\tError" > ${subdir}/import.mutational.tsv
					fi
				else
					echo "cockle-doodle-dooooooooooooooooo"
			fi
		else
			echo "HISSSSSSSSSSSSSSSSSSS"
		fi
		if [ -d ${subdir}/mutational/mutational_23S/ ]; then
			cat ${subdir}/mutational/mutational_23S/report.tsv | awk -F "\t" 'BEGIN {OFS="\t"}{if ($18 == "1") {print "sample_name",$1"("$21")"}}' | sed -e "s/sample_name/$sample_name/" > ${subdir}/import.23Sdetection.tsv
			cat ${subdir}/mutational/mutational_23S/report.tsv | awk -F "\t" 'BEGIN {OFS="\t"}{if ($18 == "1") {print "sample_name",$1"("$19")"}}' | sed -e "s/sample_name/$sample_name/" >> ${subdir}/${sample_name}.summary.23Sdetection.tsv
			if [ -s ${subdir}/${sample_name}.summary.23Sdetection.tsv ]; then
					cat ${subdir}/${sample_name}.summary.23Sdetection.tsv | echo -e "${sample_name}\t""$(cut -f 2 | paste -d, -s)" > ${subdir}/all.summary.23Sdetection.tsv
				else
					echo "heE-hAw"
			fi
			if [ ! -s ${subdir}/${sample_name}.summary.23Sdetection.tsv ]; then
				if [ -f ${subdir}/mutational/mutational_23S/log.clusters.gz ]; then
						echo -e "${sample_name}\tNone" > ${subdir}/all.summary.23Sdetection.tsv
						echo -e "${sample_name}\tNone" > ${subdir}/import.23Sdetection.tsv
					else
						echo -e "${sample_name}\tError" > ${subdir}/all.summary.23Sdetection.tsv  
						echo -e "${sample_name}\tError" > ${subdir}/import.23Sdetection.tsv
					fi
				else
					echo "cockle-doodle-dooooooooooooooooo"
			fi
		else
			echo "neigh"
		fi
		if [ -d ${subdir}/mutational/mutational_gyrA/ ]; then
			cat ${subdir}/mutational/mutational_gyrA/report.tsv | awk -F "\t" 'BEGIN {OFS="\t"}{if ($18 == "1") {print "sample_name",$1"("$19")"}}' | sed -e "s/sample_name/$sample_name/" | sed -e 's/(./(/' | sed -e 's/.)/)/' > ${subdir}/import.gyrAdetection.tsv
			cat ${subdir}/mutational/mutational_gyrA/report.tsv | awk -F "\t" 'BEGIN {OFS="\t"}{if ($18 == "1") {print "sample_name",$1"("$19")"}}' | sed -e "s/sample_name/$sample_name/" >> ${subdir}/${sample_name}.summary.gyrAdetection.tsv
			if [ -s ${subdir}/${sample_name}.summary.gyrAdetection.tsv ]; then
					cat ${subdir}/${sample_name}.summary.gyrAdetection.tsv | echo -e "${sample_name}\t""$(cut -f 2 | paste -d, -s)" > ${subdir}/all.summary.gyrAdetection.tsv
				else
					echo "heE-hAw"
			fi
			if [ ! -s ${subdir}/${sample_name}.summary.gyrAdetection.tsv ]; then
				if [ -f ${subdir}/mutational/mutational_gyrA/log.clusters.gz ]; then
						echo -e "${sample_name}\tNone" > ${subdir}/all.summary.gyrAdetection.tsv
						echo -e "${sample_name}\tNone" > ${subdir}/import.gyrAdetection.tsv
					else
						echo -e "${sample_name}\tError" > ${subdir}/all.summary.gyrAdetection.tsv  
						echo -e "${sample_name}\tError" > ${subdir}/import.gyrAdetection.tsv
					fi
				else
					echo "cockle-doodle-dooooooooooooooooo"
			fi
		else
			echo "bah"
		fi
	else
		echo "woofwoof"
	fi
	#Sort import.tsv's and join to create final import sheets for plasmidfinder, resfinder, and pointfinder/mutations
	if [ -e ${subdir}/staramr_summary.tsv ]; then
		sort -k 1 -o ${subdir}/staramr_summary.tsv ${subdir}/staramr_summary.tsv
	fi
	if [ -e ${subdir}/import.resfinder.tsv ]; then
		sort -k 1 -o ${subdir}/import.resfinder.tsv ${subdir}/import.resfinder.tsv
	fi
	if [ -e ${subdir}/import.plasmidfinder.tsv ]; then
		sort -k 1 -o ${subdir}/import.plasmidfinder.tsv ${subdir}/import.plasmidfinder.tsv
	fi
	if [ -e ${subdir}/import.mutational.tsv ]; then
		sort -k 1 -o ${subdir}/import.mutational.tsv ${subdir}/import.mutational.tsv
	fi
	if [ -e ${subdir}/all.summary.mutational.tsv ]; then
		sort -k 1 -o ${subdir}/all.summary.mutational.tsv ${subdir}/all.summary.mutational.tsv
	fi
	if [ -e ${subdir}/import.23Sdetection.tsv ]; then
		sort -k 1 -o ${subdir}/import.23Sdetection.tsv ${subdir}/import.23Sdetection.tsv
	fi
	if [ -e ${subdir}/all.summary.23Sdetection.tsv ]; then
		sort -k 1 -o ${subdir}/all.summary.23Sdetection.tsv ${subdir}/all.summary.23Sdetection.tsv
	fi
	if [ -e ${subdir}/import.gyrAdetection.tsv ]; then
		sort -k 1 -o ${subdir}/import.gyrAdetection.tsv ${subdir}/import.gyrAdetection.tsv
	fi
	if [ -e ${subdir}/all.summary.gyrAdetection.tsv ]; then
		sort -k 1 -o ${subdir}/all.summary.gyrAdetection.tsv ${subdir}/all.summary.gyrAdetection.tsv
	fi
	#Create final.summary.tsv based on organism run conditions.
	#1) Check if 23S, gyrA exist AND pointfinder does NOT exist: campy.coli, campy.upsaliensis
	mkdir -p ${subdir}import_sheets/
	if [ -e ${subdir}/import.23Sdetection.tsv ] && [ -e ${subdir}/import.gyrAdetection.tsv ] && [ ! -e ${subdir}/staramr.pointfinder.tsv ]; then
		cat ${subdir}/import.23Sdetection.tsv ${subdir}/import.gyrAdetection.tsv > ${subdir}/import.pointfinder.tsv
		join ${subdir}/staramr_summary.tsv ${subdir}/all.summary.23Sdetection.tsv > ${subdir}/determinants.tsv
		join ${subdir}/determinants.tsv ${subdir}/all.summary.23Sdetection.tsv > ${subdir}/determinants2.tsv
		join ${subdir}/determinants2.tsv ${subdir}/import.plasmidfinder.tsv | sed -e 's/\ /\t/g' > ${subdir}/final.summary.tsv
		awk '{FS=OFS="\t"}{print $0,"NARMSWF4.0"}' ${subdir}/final.summary.tsv > ${subdir}/import_sheets/final.summary."$time_stamp".tsv
		awk '{FS=OFS="\t"}{print $0,"NARMSWF4.0"}' ${subdir}/import.pointfinder.tsv > ${subdir}/import_sheets/import.pointfinder."$time_stamp".tsv
		awk '{FS=OFS="\t"}{print $0,"NARMSWF4.0"}' ${subdir}/import.plasmidfinder.tsv > ${subdir}/import_sheets/import.plasmidfinder."$time_stamp".tsv
		awk '{FS=OFS="\t"}{print $0,"NARMSWF4.0"}' ${subdir}/import.resfinder.tsv > ${subdir}/import_sheets/import.resfinder."$time_stamp".tsv
	#2) Check if 23S exists AND pointfinder exists: campy.jejuni
	elif [ -e ${subdir}/import.23Sdetection.tsv ] && [ -e ${subdir}/staramr.pointfinder.tsv ]; then
		cat ${subdir}/import.23Sdetection.tsv ${subdir}/staramr.pointfinder.tsv > ${subdir}/import.pointfinder.tsv
		join ${subdir}/staramr_summary.tsv ${subdir}/all.summary.23Sdetection.tsv > ${subdir}/determinants.tsv
		join ${subdir}/determinants.tsv ${subdir}/import.plasmidfinder.tsv | sed -e 's/\ /\t/g' > ${subdir}/final.summary.tsv
		awk '{FS=OFS="\t"}{print $0,"NARMSWF4.0"}' ${subdir}/final.summary.tsv > ${subdir}/import_sheets/final.summary."$time_stamp".tsv
		awk '{FS=OFS="\t"}{print $0,"NARMSWF4.0"}' ${subdir}/import.pointfinder.tsv > ${subdir}/import_sheets/import.pointfinder."$time_stamp".tsv
		awk '{FS=OFS="\t"}{print $0,"NARMSWF4.0"}' ${subdir}/import.plasmidfinder.tsv > ${subdir}/import_sheets/import.plasmidfinder."$time_stamp".tsv
		awk '{FS=OFS="\t"}{print $0,"NARMSWF4.0"}' ${subdir}/import.resfinder.tsv > ${subdir}/import_sheets/import.resfinder."$time_stamp".tsv
	#3) Check if mutational exists AND pointfinder does NOT exist: e.coli, vibrioc, vibriop, yersiniae
	elif [ -e ${subdir}/import.mutational.tsv ] && [ ! -e ${subdir}/staramr.pointfinder.tsv ]; then
		join ${subdir}/staramr_summary.tsv ${subdir}/all.summary.mutational.tsv > ${subdir}/determinants.tsv
		join ${subdir}/determinants.tsv ${subdir}/import.plasmidfinder.tsv | sed -e 's/\ /\t/g' > ${subdir}/final.summary.tsv
		awk '{FS=OFS="\t"}{print $0,"NARMSWF4.0"}' ${subdir}/final.summary.tsv > ${subdir}/import_sheets/final.summary."$time_stamp".tsv
		awk '{FS=OFS="\t"}{print $0,"NARMSWF4.0"}' ${subdir}/import.mutational.tsv > ${subdir}/import_sheets/import.pointfinder."$time_stamp".tsv
		awk '{FS=OFS="\t"}{print $0,"NARMSWF4.0"}' ${subdir}/import.plasmidfinder.tsv > ${subdir}/import_sheets/import.plasmidfinder."$time_stamp".tsv
		awk '{FS=OFS="\t"}{print $0,"NARMSWF4.0"}' ${subdir}/import.resfinder.tsv > ${subdir}/import_sheets/import.resfinder."$time_stamp".tsv
	#4) Check if mutational does NOT exist AND pointfinder does exist: Salmonella
	elif [ ! -e ${subdir}/import.mutational.tsv ] && [ -e ${subdir}/staramr.pointfinder.tsv ]; then
		join ${subdir}/staramr_summary.tsv ${subdir}/import.plasmidfinder.tsv | sed -e 's/\ /\t/g' > ${subdir}/final.summary.tsv
		awk '{FS=OFS="\t"}{print $0,"NARMSWF4.0"}' ${subdir}/final.summary.tsv > ${subdir}/import_sheets/final.summary."$time_stamp".tsv
		awk '{FS=OFS="\t"}{print $0,"NARMSWF4.0"}' ${subdir}/staramr.pointfinder.tsv > ${subdir}/import_sheets/import.pointfinder."$time_stamp".tsv
		awk '{FS=OFS="\t"}{print $0,"NARMSWF4.0"}' ${subdir}/import.plasmidfinder.tsv > ${subdir}/import_sheets/import.plasmidfinder."$time_stamp".tsv
		awk '{FS=OFS="\t"}{print $0,"NARMSWF4.0"}' ${subdir}/import.resfinder.tsv > ${subdir}/import_sheets/import.resfinder."$time_stamp".tsv
	#5) Check if mutational does NOT exist and pointfinder does NOT exist: Other
	elif [ ! -e ${subdir}/import.mutational.tsv ] && [ ! -e ${subdir}/staramr.pointfinder.tsv ]; then
		join ${subdir}/staramr_summary.tsv ${subdir}/import.plasmidfinder.tsv | sed -e 's/\ /\t/g' > ${subdir}/final.summary.tsv
		awk '{FS=OFS="\t"}{print $0,"NARMSWF4.0"}' ${subdir}/final.summary.tsv > ${subdir}/import_sheets/final.summary."$time_stamp".tsv
		awk '{FS=OFS="\t"}{print $0,"NARMSWF4.0"}' ${subdir}/import.plasmidfinder.tsv > ${subdir}/import_sheets/import.plasmidfinder."$time_stamp".tsv
		awk '{FS=OFS="\t"}{print $0,"NARMSWF4.0"}' ${subdir}/import.resfinder.tsv > ${subdir}/import_sheets/import.resfinder."$time_stamp".tsv
	else
		"idk man didn't work. message Justin to see if he can troubleshoot wtf happened. this script is way too long."
	fi
done
