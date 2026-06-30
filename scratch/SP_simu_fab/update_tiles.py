import os
import glob
import re

def load_numbers_from_file(file_path):
    """Loads the XXX-XXX numbers from the input text file."""
    with open(file_path, 'r') as f:
        return set(line.strip() for line in f if re.match(r'\d{3}-\d{3}', line.strip()))

def find_matching_directories(base_path):
    """Finds all directories matching the required pattern using wildcards."""
    search_pattern = os.path.join(base_path, "output_*/run_sp_tile_ngmix_Ng1u_*/ngmix_runner/output")
    return glob.glob(search_pattern)

def extract_numbers_from_filenames(directory):
    """Extracts XXX-XXX patterns from filenames in the given directory."""
    pattern = re.compile(r'ngmix-(\d{3}-\d{3})\.fits')
    found_numbers = set()
    
    for filename in os.listdir(directory):
        match = pattern.search(filename)
        if match:
            found_numbers.add(match.group(1))
    
    return found_numbers

def write_output(file_path, data):
    """Writes a set of numbers to a file."""
    with open(file_path, 'w') as f:
        for item in sorted(data):
            f.write(item + '\n')

def main(input_file, base_directory, found_output, not_found_output):
    """Main function to process the numbers and check their existence in files."""
    input_numbers = load_numbers_from_file(input_file)
    matching_dirs = find_matching_directories(base_directory)
    
    available_numbers = set()
    for directory in matching_dirs:
        available_numbers.update(extract_numbers_from_filenames(directory))
    
    found = input_numbers & available_numbers
    not_found = input_numbers - available_numbers
    
    write_output(found_output, found)
    write_output(not_found_output, not_found)
    
    print(f"Processing complete. Found: {len(found)}, Not Found: {len(not_found)}")

if __name__ == "__main__":
    input_file = "/home/hervas/n25/SP_simu_fab/tile_numbers_real.txt"  # Change this to your actual input file path
    base_directory = "./outputs"  # Base directory to start searching
    found_output = "found.txt"
    not_found_output = "numbers_run.txt"
    
    main(input_file, base_directory, found_output, not_found_output)
