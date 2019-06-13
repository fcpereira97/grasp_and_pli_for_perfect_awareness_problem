#!/bin/bash

for filename in instances/groupA/*; do
	echo "$filename"
	g++ -o grasp-pa grasp-pa.cpp && ./grasp-pa 2 1 2 1 1800 "$filename" output/output_random_plus_greedy_groupA.csv >> logs/log_random_plus_greedy_groupA.txt
done