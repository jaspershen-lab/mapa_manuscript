#!/usr/bin/env bash

# The file containing the URLs
LINKS_FILE="proteomics_links.txt"

# The directory where you want to save the downloaded files
DOWNLOAD_DIR="proteomics"

# Create the download directory if it doesn't exist
mkdir -p "$DOWNLOAD_DIR"

# Read file line by line, ensuring the last line is processed even if it lacks a trailing newline
while IFS= read -r line || [ -n "$line" ]; do
	# Remove any trailing newline or carriage return characters
	url="${line//[$'\r\n']}"
		    
	# Skip if the line is empty (after trimming)
	[[ -z "$url" ]] && continue
			    
	echo "Downloading $url into $DOWNLOAD_DIR ..."
	wget -c -P "$DOWNLOAD_DIR" "$url"

done < "$LINKS_FILE"
