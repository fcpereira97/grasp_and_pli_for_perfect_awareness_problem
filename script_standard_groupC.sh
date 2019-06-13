#!/bin/bash

for filename in instances/groupC/*; do
	echo "$filename"
	g++ -o grasp-pa grasp-pa.cpp && ./grasp-pa 1 2 3600 "$filename" output/output_standard_groupC.csv >> logs/log_standard_groupC.txt
done