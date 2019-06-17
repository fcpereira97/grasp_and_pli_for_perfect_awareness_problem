#!/bin/bash

for filename in instances/groupC/*; do
	echo "$filename"
	g++ -o grasp-pa grasp-pa.cpp && ./grasp-pa 3 100 2 0.05 3600 "$filename" output/output_sampled_greedy_groupC.csv >> logs/log_sampled_greedy_groupC.txt
done