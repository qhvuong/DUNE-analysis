import os

# Generate the expected filenames
expected_files = [f"MiniRun3_1E19_RHC.flow_v6.{d:05d}.FLOW.h5" for d in range(1024)]


# Get a list of all files in the directory
files_in_directory = os.listdir("/pnfs/dune/tape_backed/users/mkramer/prod/MiniRun3/MiniRun3_1E19_RHC/MiniRun3_1E19_RHC.flow_v6/")

# Find the missing files
missing_files = [filename for filename in expected_files if filename not in files_in_directory]

# Print the missing files
print("Missing files:")
for filename in missing_files:
    print(filename)

