import os
import argparse
import subprocess
import glob
import pandas as pd
import matplotlib.pyplot as plt

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

#function to graph the results into R graph
def plotfunc(filename, output_dir):
    df = pd.read_csv(filename, sep='\t', header=None, usecols=[0,1], names=['Chromosome', ['Position']])
    df['Chromosome'] = df['Chromosome'].astype(str) #converts chromosome into strings for categorical data (i.e. chromosome 1 will be 1)

    #columns chrom and pos
    plot = df.plot(kind = 'scatter', x = 'Chromosome', y = 'Position', title=filename)
    
    #save plot to output directory
    prefix_plot = os.path.join(output_dir, os.path.basename(filename).replace(".tab", "_plot.jpg"))
    plt.savefig(prefix_plot)
    print(f"Plot saved as {prefix_plot}")
    plt.close()


def main(input_directory, output_dir):
    #debug output directory
    print(f"Output directory: {output_dir}")
    #create output directory if directory does not exist
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    #create annotated files
    script = "/mnt/storage8/mtan/scripts/annotate.sh" #to change where bash script annotate.sh is stored
    input_path = [os.path.join(input_directory)]
    os.chdir(input_path) #changing into the directory where the vcf files are there

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
            plot = plotfunc(os.path.join(annotation_path, tab_file))
    
    



if __name__=="__main__":
    # Set up argparse to handle directory input and optional comparison
    parser = argparse.ArgumentParser(description="Going from annotated Pf genomes to plot graph showing SNP locations in chromosomes.")
    parser.add_argument("input_directory", help="Directory containing the vcf files.")
    parser.add_argument("--output", default=".", help="Directory where output plots will be saved (default: current directory)")
    
    args = parser.parse_args()

    main(args.input_directory, args.output)

