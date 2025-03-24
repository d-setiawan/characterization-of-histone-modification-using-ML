#!/bin/bash

# Loop through each gzip file and unzip it
for file in *.gz; do
    echo "Unzipping $file..."
    gunzip "$file"
done

echo "All files unzipped!"