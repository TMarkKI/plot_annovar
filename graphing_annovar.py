import os
import argparse
import subprocess
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Try importing r, install it if not found
#try:
#    import rpy2
#    import rpy2.objects as robjects
#except ImportError:
#    subprocess.check_call(["pip", "install", "rpy2"])
#    import rpy2
#    import rpy2.objects as robjects

#taken from https://stackoverflow.com/questions/2651874/embed-bash-in-python by Ian Bicking to run the annotate.sh bash script
def run_script_annovar(script, stdin=None):
    process = subprocess.Popen(['bash', '-c', script],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        stdin=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode:
            raise ScriptException(process.returncode, stdout, stderr, script)
    return stdout, stderr

class ScriptException(Exception):
      def __init__(self, returncode, stdout, stderr, script):
            self.returncode = returncode
            self.stdout = stdout
            self.stderr = stderr
            Exception().__init__('Error in annovar script')

#process annotated table to seek out discrepancies in the samples; uses positional/column options so this will affect the columns specific to the number of samples
def discrepancy_finder(table):
    df = pd.read_csv(table, sep='\t', header=None)
     
    #conditional insertion for comparison between the genotypes of two samples
    df['16h_genotype_score'] = np.where(df[5] == '{0/0}', 0,
                                        np.where(df[5] == '{0/1}', 1.5, 
                                                 np.where(df[5] == '{1/1}', 1, np.nan)))

    df['8h_genotype_score'] = np.where(df[11] == '{0/0}', 0,
                                       np.where(df[11] == '{0/1}', 1.5, 
                                                np.where(df[11] == '{1/1}', 1, np.nan)))

    df['discrepant score'] = df['16h genotype score'] + df['8h genotype score'] #need to tabulate score to assign condition of gt

    df['condition'] = np.where(df['discrepant_score'] == 0, 'wild-type',
                               np.where(df['discrepant_score'] == 1.5, 'wild-type and mixed',
                                        np.where(df['discrepant_score'] == 3, 'mixed',
                                                 np.where(df['discrepant_score'] == 2.5, 'alternate and mixed',
                                                          np.where(df['discrepant_score'] == 2, 'alternate', 
                                                                   'alternate and wild-type')))))
    return df



#function to graph the results into R graph
def plotfunc(filename, output_dir):
    df = pd.read_csv(filename, sep='\t')
    df_plot = discrepancy_finder(df)
    #df = pd.read_csv(filename, sep='\t', header=None, usecols=[0,1,15], names=['Chromosome', 'Position','condition'])
    #df['Chromosome'] = df['Chromosome'].astype(str) #converts chromosome into strings for categorical data (i.e. chromosome 1 will be 1)

    #columns chrom and pos; added another section to determine colour based on discrepancy/acceptancy
    #colour conditions
    conditions_colour_map = {
        'wild-type': 'blue',
        'wild-type and mixed': 'green',
        'mixed': 'yellow',
        'alternate and mixed': 'orange',
        'alternate': 'red',
        'alternate and wild-type': 'purple'
    }

    df_plot['colour'] = df_plot['condition'].map(conditions_colour_map)
    plot = df_plot.plot(kind='scatter', x='Chromosome', y='Position', color=df_plot['colour'], title=filename)


    #df_plot = pd.read_csv(filename, sep='\t', header=None, usecols=[0,1,15], names=['Chromosome', 'Position', 'condition'])
    #plot = df_plot.plot(kind = 'scatter', x = 'Chromosome', y = 'Position', title=filename)
    
    #save plot to output directory
    prefix_plot = os.path.join(output_dir, os.path.basename(filename).replace(".tab", "_plot.jpg"))
    plt.savefig(prefix_plot)
    print(f"Plot saved as {prefix_plot}")
    plt.close()



#main section of code
def main(input_directory, output_dir):
    #debug output directory
    print(f"Output directory: {output_dir}")
    #create output directory if directory does not exist
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    #create annotated files
    script = "/mnt/storage8/mtan/scripts/annotate.sh" #to change where bash script annotate.sh is stored
    #input_path = [os.path.join(input_directory)]
    os.chdir(input_directory) #changing into the directory where the vcf files are there

    try:
        stdout, stderr = run_script_annovar(script)
        print("Standard Output:", stdout.decode())
        print("Standard Error:", stderr.decode())
    except ScriptException as e:
        print(f"Script failed with return code {e.returncode}")
        print(f"Error output: {e.stderr.decode()}")
    
    #changing the directory to the annotate directory that the annotated vcf and tab files are stored in
    annotation_path = os.path.join(input_directory, "annotate")
    os.chdir(annotation_path)
    
    tab_directory = glob.glob("*.tab")
    if not tab_directory:
         print("No .tab files, recheck directory")
    else:
         for tab_file in tab_directory:
            print(f"Processing {tab_file}")
            plotfunc(os.path.join(annotation_path, tab_file), output_dir)
    
    



if __name__=="__main__":
    # Set up argparse to handle directory input and optional comparison
    parser = argparse.ArgumentParser(description="Going from annotated Pf genomes to plot graph showing SNP locations in chromosomes.")
    parser.add_argument("input_directory", help="Directory containing the vcf files.")
    parser.add_argument("--output", default=".", help="Directory where output plots will be saved (default: current directory)")
    
    args = parser.parse_args()

    main(args.input_directory, args.output)

