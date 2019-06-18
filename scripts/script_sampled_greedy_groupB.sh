#!/bin/bash

for filename in instances/groupB/*; do
	echo "$filename"
	g++ -o grasp-pa grasp-pa.cpp && ./grasp-pa 3 1 2 1 600 "$filename" output/output_sampled_greedy_groupB.csv >> logs/log_sampled_greedy_groupB.txt
done