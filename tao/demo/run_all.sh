#!/bin/bash

# 1. Define your tasks: "Folder Name | Command"
# Keep the pipe symbol (|) as a separator


tao -Q -V -N31 -s1e12 -B1 -D0/300e3 -T20e3
tao -Qexample_B_load -N51 -B2 -p1e12 -D0/1000e3
TASKS=(
    "tao example_C"
    "tao example_D -P2"
    "tao example_E"
    "example_F.job"
    "tao example_G"
)


echo "Launching unique parallel tasks..."

# 2. Loop through the array
for task in "${TASKS[@]}"; do
            echo "--- Running '$task'"
            eval "$task" 
            echo "--- $task Done! ---"
done

# 3. Wait for all background tasks
jobs
wait

echo "---------------------------------------"
echo "ALL TASKS FINISHED."
echo 

