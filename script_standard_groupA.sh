#!/bin/bash

for filename in instances/groupA/*; do
	echo "$filename"
	g++ -o grasp-pa grasp-pa.cpp && ./grasp-pa 1 2 1800 "$filename" output/output_standard_groupA.csv >> logs/log_standard_groupA.txt
done