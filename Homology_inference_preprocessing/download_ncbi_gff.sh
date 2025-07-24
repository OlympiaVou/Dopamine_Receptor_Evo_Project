echo " ======================================================"
echo "|                                                      |"
echo " ~~~~~~~~~~~~~~~~~~~~~~~~START~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "|                                                      |"
echo -e " ======================================================\n"


GCF_IDs=$(python3 find_GCF_column.py ${1})

for proteome in ${GCF_IDs}
do
	var=$(echo "${proteome}" |grep -oP "(?<=GCF_)\d+")
	first3="${var:0:3}"
	second3="${var:3:3}"
	last3="${var:6:3}"

	# Save the link of the GCF to search for the gff file
	base_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${first3}/${second3}/${last3}/"
	echo ${base_url}
        echo -e "------------------------------------------------------"
	echo -e "Searching for GFF for ${proteome}...\n"

	# Get the HTML page for the subdirectory
	page=$(curl -s "${base_url}/")

	# Extract the full subfolder name
	folder=$(echo "$page" | grep -oP "GCF_${var}\.\d+_[^/]+" | head -1)

	# Check if the subfolder exists
	if [ -z "$folder" ]; then
		echo "Could not find subfolder for ${proteome}"
		continue
	fi

	# Save the complete link to download the gff
	gff_url="${base_url}/${folder}/${folder}_genomic.gff.gz"
	echo -e "Downloading: $gff_url\n"

        # Download the gff file
	wget "$gff_url"

	echo -e "------------------------------------------------------\n"
done

echo "======================================================"
echo "|                                                    |"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~END~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "|                                                    |"
echo "======================================================"

