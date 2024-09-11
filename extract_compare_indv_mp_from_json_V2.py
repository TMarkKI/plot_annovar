#version 2 now contains an output option and for the counts to be saved in a txt

import os
import json
import argparse
import pandas as pd

# Try importing alive-progress, install it if not found
try:
    from alive_progress import alive_bar
except ImportError:
    import subprocess
    subprocess.check_call(["pip", "install", "alive-progress"])
    from alive_progress import alive_bar

def process_json_file(json_file, output_dir):
    # Read and parse the JSON file
    with open(json_file, 'r') as file:
        json_data = json.load(file)

    print(f"Processing {json_file}:")

    # Process the variants and add a source column
    DR_vars = pd.DataFrame(json_data.get('dr_variants', []))
    DR_vars['source'] = 'DR'

    other_vars = pd.DataFrame(json_data.get('other_variants', []))
    other_vars['source'] = 'other'

    fail_vars = pd.DataFrame(json_data.get('fail_variants', []))
    fail_vars['source'] = 'fail'

    # Save the dataframes to CSV files in the specified output directory
    output_prefix = os.path.basename(json_file).replace('.fastq.gz.results.json', '')
    DR_vars.to_csv(os.path.join(output_dir, f"DR_vars_{output_prefix}.csv"), index=False)
    other_vars.to_csv(os.path.join(output_dir, f"other_vars_{output_prefix}.csv"), index=False)
    fail_vars.to_csv(os.path.join(output_dir, f"fail_vars_{output_prefix}.csv"), index=False)

    # Combine all variants into one DataFrame with the source column
    compiled_vars = pd.concat([DR_vars, other_vars, fail_vars], ignore_index=True)
    compiled_vars['sample'] = output_prefix
    compiled_vars.to_csv(os.path.join(output_dir, f"compiled_{output_prefix}.csv"), index=False)
    print(compiled_vars)
    return compiled_vars

def compare_variants(sample1_df, sample2_df, columns_to_compare=['chrom', 'pos', 'ref', 'alt']):
    # Prepare a list to store comparison results
    comparison_results = []

    # Iterate through each variant in sample1_df
    for idx1, row1 in sample1_df.iterrows():
        found_in_sample2 = False
        for idx2, row2 in sample2_df.iterrows():
            # Check if only the specified columns match between the two rows
            if row1[columns_to_compare].equals(row2[columns_to_compare]):
                found_in_sample2 = True
                if row1['source'] == row2['source']:
                    comparison_results.append((row1.to_dict(), 'Y'))  # Same variant and same status
                else:
                    comparison_results.append((row1.to_dict(), '?'))  # Same variant but different status
                break
        
        if not found_in_sample2:
            comparison_results.append((row1.to_dict(), 'N'))  # Variant not found in sample2

    # Check for variants in sample2 that are not in sample1
    for idx2, row2 in sample2_df.iterrows():
        found_in_sample1 = False
        for idx1, row1 in sample1_df.iterrows():
            if row2[columns_to_compare].equals(row1[columns_to_compare]):
                found_in_sample1 = True
                if row2['source'] == row1['source']:
                    comparison_results.append((row2.to_dict(), 'Y'))
                else:
                    comparison_results.append((row2.to_dict(), '?'))
                break
            
        if not found_in_sample1:
            comparison_results.append((row2.to_dict(), 'N'))  # Variant not found in sample1

    # Convert the results to a DataFrame
    comparison_df = pd.DataFrame([dict(**comp[0], comparison=comp[1]) for comp in comparison_results])

    return comparison_df

def main(directory, output_dir, compare_csv=None):
    # Debug: Print the output directory
    print(f"Output directory: {output_dir}")

    # Create the output directory if it does not exist
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # List all JSON files matching the pattern in the specified directory
    json_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".fastq.gz.results.json")]

    # Process each JSON file
    processed_dataframes = {}
    with alive_bar(len(json_files), title="Processing JSON Files") as bar:
        for json_file in json_files:
            processed_dataframes[json_file] = process_json_file(json_file, output_dir)
            bar()

    # If comparison is enabled, perform pairwise comparisons
    if compare_csv:
        comparison_pairs = pd.read_csv(compare_csv)
        with alive_bar(len(comparison_pairs), title="Performing comparisons") as bar:
            for _, row in comparison_pairs.iterrows():
                sample1_file = 'compiled_' + row['sample1'] + '.csv'
                sample2_file = 'compiled_' + row['sample2'] + '.csv'

                # Look for sample1 and sample2 compiled CSV files in the output directory
                sample1_df = pd.read_csv(os.path.join(output_dir, sample1_file))
                sample2_df = pd.read_csv(os.path.join(output_dir, sample2_file))

                if sample1_df is not None and sample2_df is not None:
                    comparison_result = compare_variants(sample1_df, sample2_df)

                    # Count Y, ?, and N in the comparison results of sample1_df
                    y_count1 = ((comparison_result['comparison'] == 'Y') & (comparison_result['sample'] == row['sample1'])).sum()
                    q_count1 = ((comparison_result['comparison'] == '?') & (comparison_result['sample'] == row['sample1'])).sum()
                    n_count1 = ((comparison_result['comparison'] == 'N') & (comparison_result['sample'] == row['sample1'])).sum()
                    y_DR_count1 = ((comparison_result['comparison'] == 'Y') & (comparison_result['source'] == 'DR') & (comparison_result['sample'] == row['sample1'])).sum()
                    q_DR_count1 = ((comparison_result['comparison'] == '?') & (comparison_result['source'] == 'DR') & (comparison_result['sample'] == row['sample1'])).sum()
                    n_DR_count1 = ((comparison_result['comparison'] == 'N') & (comparison_result['source'] == 'DR') & (comparison_result['sample'] == row['sample1'])).sum()
                    y_other_count1 = ((comparison_result['comparison'] == 'Y') & (comparison_result['source'] == 'other') & (comparison_result['sample'] == row['sample1'])).sum()
                    q_other_count1 = ((comparison_result['comparison'] == '?') & (comparison_result['source'] == 'other') & (comparison_result['sample'] == row['sample1'])).sum()
                    n_other_count1 = ((comparison_result['comparison'] == 'N') & (comparison_result['source'] == 'other') & (comparison_result['sample'] == row['sample1'])).sum()

                    # Count Y, ?, and N in the comparison results of sample2_df
                    y_count2 = ((comparison_result['comparison'] == 'Y') & (comparison_result['sample'] == row['sample2'])).sum()
                    q_count2 = ((comparison_result['comparison'] == '?') & (comparison_result['sample'] == row['sample2'])).sum()
                    n_count2 = ((comparison_result['comparison'] == 'N') & (comparison_result['sample'] == row['sample2'])).sum()
                    y_DR_count2 = ((comparison_result['comparison'] == 'Y') & (comparison_result['source'] == 'DR') & (comparison_result['sample'] == row['sample2'])).sum()
                    q_DR_count2 = ((comparison_result['comparison'] == '?') & (comparison_result['source'] == 'DR') & (comparison_result['sample'] == row['sample2'])).sum()
                    n_DR_count2 = ((comparison_result['comparison'] == 'N') & (comparison_result['source'] == 'DR') & (comparison_result['sample'] == row['sample2'])).sum()
                    y_other_count2 = ((comparison_result['comparison'] == 'Y') & (comparison_result['source'] == 'other') & (comparison_result['sample'] == row['sample2'])).sum()
                    q_other_count2 = ((comparison_result['comparison'] == '?') & (comparison_result['source'] == 'other') & (comparison_result['sample'] == row['sample2'])).sum()
                    n_other_count2 = ((comparison_result['comparison'] == 'N') & (comparison_result['source'] == 'other') & (comparison_result['sample'] == row['sample2'])).sum()

                    # Print the summary to the command line
                    #print(f"Comparison between {os.path.basename(sample1_file)} and {os.path.basename(sample2_file)}:")
                    #print(f"For {row['sample1']} Y: {y_count1}, ?: {q_count1}, N: {n_count1}")
                    #print(f"For {row['sample1']} DR - Y: {y_DR_count1}, ?: {q_DR_count1}, N: {n_DR_count1}")
                    #print(f"For {row['sample1']} Other - Y: {y_other_count1}, ?: {q_other_count1}, N: {n_other_count1}")

                    #print(f"For {row['sample2']} Y: {y_count2}, ?: {q_count2}, N: {n_count2}")
                    #print(f"For {row['sample2']} DR - Y: {y_DR_count2}, ?: {q_DR_count2}, N: {n_DR_count2}")
                    #print(f"For {row['sample2']} Other - Y: {y_other_count2}, ?: {q_other_count2}, N: {n_other_count2}")


                    # Summary output for txt file
                    summary = []
                    summary.append(f"Comparison between {os.path.basename(sample1_file)} and {os.path.basename(sample2_file)}:")
                    summary.append(f"For {row['sample1']} Y: {y_count1}, ?: {q_count1}, N: {n_count1}")
                    summary.append(f"For {row['sample1']} DR - Y: {y_DR_count1}, ?: {q_DR_count1}, N: {n_DR_count1}")
                    summary.append(f"For {row['sample1']} Other - Y: {y_other_count1}, ?: {q_other_count1}, N: {n_other_count1}")
                    summary.append(f"For {row['sample2']} Y: {y_count2}, ?: {q_count2}, N: {n_count2}")
                    summary.append(f"For {row['sample2']} DR - Y: {y_DR_count2}, ?: {q_DR_count2}, N: {n_DR_count2}")
                    summary.append(f"For {row['sample2']} Other - Y: {y_other_count2}, ?: {q_other_count2}, N: {n_other_count2}")
                    summary.append("\n")

                    #Display summary on command line
                    print("\n".join(summary))

                    # Save the comparison result as CSV in the output directory
                    comparison_output_prefix = f"comparison_{row['sample1']}_vs_{row['sample2']}"
                    with open(os.path.join(output_dir, f"{comparison_output_prefix}.txt"), 'w') as f:
                        f.write("\n".join(summary)) #txt
                    comparison_result.to_csv(os.path.join(output_dir, f"{comparison_output_prefix}.csv"), index=False) #csv
                else:
                    print(f"Error: One or both files not found for comparison: {row['sample1']} vs {row['sample2']}")
                bar()

if __name__ == "__main__":
    # Set up argparse to handle directory input and optional comparison
    parser = argparse.ArgumentParser(description="Process JSON files in a directory with optional pairwise comparison.")
    parser.add_argument("directory", help="Directory containing the JSON files")
    parser.add_argument("--output", default=".", help="Directory where output CSV files will be saved (default: current directory)")
    parser.add_argument("--compare", help="CSV file specifying pairwise comparisons (two columns: sample1, sample2)")

    args = parser.parse_args()

    # Run the main function with the provided directory, output directory, and comparison CSV (if provided)
    main(args.directory, args.output, args.compare)