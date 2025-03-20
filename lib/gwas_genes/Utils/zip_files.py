import os
import shutil
import zipfile

def copy_and_zip_csvs(csv_files, output_dir, zip_filename="csv_files.zip"):
    """
    Copies the CSV files to the output_dir, zips the directory, and returns the zip file path.

    Args:
        csv_files (list): List of CSV file paths.
        output_dir (str): Directory where the CSV files will be copied.
        zip_filename (str): Name of the resulting zip file (default: "csv_files.zip").

    Returns:
        str: Path to the created zip file.
    """
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Copy each CSV file to the output directory
    for csv_file in csv_files:
        if os.path.isfile(csv_file):
            shutil.copy(csv_file, output_dir)
        else:
            print(f"Warning: {csv_file} does not exist or is not a file.")

    # Determine the zip file path (placing it one level above the output_dir)
    parent_dir = os.path.dirname(output_dir)
    zip_filepath = os.path.join(parent_dir, zip_filename)

    # Create a zip file and add all files from the output directory
    with zipfile.ZipFile(zip_filepath, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, _, files in os.walk(output_dir):
            for file in files:
                file_path = os.path.join(root, file)
                # Use relative path for archiving
                arcname = os.path.relpath(file_path, start=parent_dir)
                zipf.write(file_path, arcname)

    return zip_filepath

if __name__ == "__main__":
    # Example list of CSV file paths (update these paths as needed)
    csv_files = ["file1.csv", "file2.csv", "file3.csv"]
    output_dir = "output_csvs"  # Directory where CSV files will be copied

    zip_path = copy_and_zip_csvs(csv_files, output_dir)
    print("Zip file created at:", zip_path)


