import argparse
import yaml
import snakemake
import os

parser = argparse.ArgumentParser()
parser.add_argument("--sample-json", type=str, required=True, help="Path to a jason file that contains sample information")
parser.add_argument("--num-jobs", type=int, default=1, help="Number of jobs to run in parallel")
parser.add_argument("--fa", type=str, required=True, help="Path to genome fasta file for bowtie2 and STAR indexing")
parser.add_argument("--gtf", type=str, required=True, help="Path to genome gtf file for STAR indexing")
parser.add_argument("--uid", type=int, required=True, help="UID on the host server, use `id -u` command to find out")
parser.add_argument("--gid", type=int, required=True, help="GID on the host server, use `id -g` command to find out")

args = parser.parse_args()

# Generate the configuration file based on the command-line arguments
CONFIG = {
    "num_jobs": args.num_jobs,
    "SAMPLES_JSON": args.sample_json,
    "fasta":args.fa,
    "gtf":args.gtf,
    "dir_names":{
      "DNA_index_dir": "outputs/0_index/bowtie2/user",
      "RNA_index_dir": "outputs/0_index/BWA/",
      "decode_dir": "outputs/1_decode/",
      "mapped_DNA_dir": "outputs/2_DNAmapping/",
      "mapped_RNA_dir": "outputs/3_RNAmapping/",
      "dedup_dir": "outputs/4_dedup/",
      "merge_dir": "outputs/5_merge/",
      "stats_dir": "outputs/6_stats/"
    }
}

#with open("config.yaml", "w") as f:
#    yaml.dump(CONFIG, f)
#
# # Run the Snakemake pipeline using the configuration file
snakemake.snakemake(snakefile="Snakefile", config=CONFIG, cores=args.num_jobs, printshellcmds=True) # 

# change output file permissions such that users from the host machine can have permissions.
folder_path = '/snakemake/outputs/'

for root, dirs, files in os.walk(folder_path):
    # Change the ownership and permissions of the folder
    os.chown(root, args.uid, args.gid)
    #os.chmod(root, 0755)

    # Change the ownership and permissions of each file
    for file in files:
        file_path = os.path.join(root, file)
        os.chown(file_path, args.uid, args.gid)
        #os.chmod(file_path, 0755)
