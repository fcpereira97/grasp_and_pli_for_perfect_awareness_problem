#!/bin/bash

for filename in instances/groupC/*; do
	echo "$filename"
	g++ -o grasp-pa grasp-pa.cpp && ./grasp-pa 1 100 2 0.05 1800 "$filename" output/output_standard_groupC.csv >> logs/log_standard_groupC.txt
done