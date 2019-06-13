#!/bin/bash

for filename in instances/groupC/*; do
	echo "$filename"
	g++ -o grasp-pa grasp-pa.cpp && ./grasp-pa 2 2 3600 "$filename" output/output_random_plus_greedy_groupC.csv >> logs/log_random_plus_greedy_groupC.txt
done