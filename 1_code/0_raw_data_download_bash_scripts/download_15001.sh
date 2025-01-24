#!/usr/bin/env bash

# The file containing the FTP URLs
LINKS_FILE="links_CRA015001_snrna_seq.txt"

# The directory where you want to save the downloaded files
DOWNLOAD_DIR="CRA015001_snrna_seq"

# Create the download directory if it doesn't exist
mkdir -p "$DOWNLOAD_DIR"

# Read file line by line
while IFS= read -r line; do
	# Remove any trailing newline or carriage return characters
	ftp_url="${line//[$'\r\n']}"
	# Skip if the line is empty
	[[ -z "$ftp_url" ]] && continue
	echo "Downloading $ftp_url into $DOWNLOAD_DIR ..."
	# Download the file into the specified folder
	wget -P "$DOWNLOAD_DIR" "$ftp_url"

done < "$LINKS_FILE"

