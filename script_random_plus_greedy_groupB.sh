#!/bin/bash

for filename in instances/groupB/*; do
	echo "$filename"
	g++ -o grasp-pa grasp-pa.cpp && ./grasp-pa 2 2 1800 "$filename" output/output_random_plus_greedy_groupB.csv >> logs/log_random_plus_greedy_groupB.txt
done