#!/bin/bash

# 1. Define your tasks: "Folder Name | Command"
# Keep the pipe symbol (|) as a separator
TASKS=(
    "3_flexure | flex"
    "4_river_transport | conic_island"
    "5_lakes | 1lake"
    "5_lakes | 2lakes"
    "6_tectonics | demo"
    "7_climate | sinusoidal"
    "8_Ebro_evol | Ebro_evol"
    "9_Iberian_drainage | Iberia_present"
)

echo "Launching unique parallel tasks..."

# 2. Loop through the array
for task in "${TASKS[@]}"; do
    # Split the string into folder and command using 'IFS' (Internal Field Separator)
    folder=$(echo "$task" | cut -d'|' -f1 | xargs)
    command=$(echo "$task" | cut -d'|' -f2- | xargs)

    if [ -d "$folder" ]; then
        (
            cd "$folder" || exit
            echo "[In $folder]: Running '$command'"
            
            eval "tisc $command" > $command.log
            
            echo "--- [In $folder]: Done! ---"
        ) &
    else
        echo "Skipping: Directory '$folder' not found."
    fi
done

# 3. Wait for all background tasks
jobs
wait

echo "---------------------------------------"
echo "ALL TASKS FINISHED."
echo 

