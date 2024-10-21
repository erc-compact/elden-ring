#!/usr/bin/env python3
import sys, os, subprocess
import argparse
import csv
import glob
import stat
import time
from multiprocessing import Pool

def find_candidate_details(csv_file, candidate):
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if candidate in row['png_path']:
                # Store the details in a dictionary
                candidate_details = {
                    'pointing_id': row['pointing_id'],
                    'beam_id': row['beam_id'],
                    'beam_name': row['beam_name'],
                    'source_name': row['source_name'],
                    'ra': row['ra'],
                    'dec': row['dec'],
                    'gl': row['gl'],
                    'gb': row['gb'],
                    'mjd_start': row['mjd_start'],
                    'utc_start': row['utc_start'],
                    'f0_user': row['f0_user'],
                    'f0_opt': row['f0_opt'],
                    'f0_opt_err': row['f0_opt_err'],
                    'f1_user': row['f1_user'],
                    'f1_opt': row['f1_opt'],
                    'f1_opt_err': row['f1_opt_err'],
                    'acc_user': row['acc_user'],
                    'acc_opt': row['acc_opt'],
                    'acc_opt_err': row['acc_opt_err'],
                    'dm_user': row['dm_user'],
                    'dm_opt': row['dm_opt'],
                    'dm_opt_err': row['dm_opt_err'],
                    'sn_fft': row['sn_fft'],
                    'sn_fold': row['sn_fold'],
                    'pepoch': row['pepoch'],
                    'maxdm_ymw16': row['maxdm_ymw16'],
                    'dist_ymw16': row['dist_ymw16'],
                    'pics_trapum_ter5': row['pics_trapum_ter5'],
                    'pics_palfa': row['pics_palfa'],
                    'pics_meerkat_l_sband_combined_best_recall': row['pics_meerkat_l_sband_combined_best_recall'],
                    'pics_palfa_meerkat_l_sband_best_fscore': row['pics_palfa_meerkat_l_sband_best_fscore'],
                    'png_path': row['png_path'],
                    'metafile_path': row['metafile_path'],
                    'filterbank_path': row['filterbank_path'],
                    'candidate_tarball_path': row['candidate_tarball_path']
                }
                return candidate_details
    return None

def submit_job(command):
    try:
        subprocess.run(command, shell=True, check=True)
        print(True, command)
    except subprocess.CalledProcessError as e:
        print(f"Error in running command: {command}")
        print(e)

def create_submit_file(work_dir, cand, csv_file, fil, num_cores):
    
    start_time = time.time()
    
    candidate_details = find_candidate_details(csv_file, cand)
    if not candidate_details:
        print(f"Candidate details not found for {cand}")
        return

    dir_name = f"{cand.split('/')[0]}_{cand.split('/')[1]}_{cand.split('/')[-1].split('.')[1].split('_')[-1]}"
    name = cand.split('/')[-1].split('.')[1].split('_')[-1]  # This is the name of the candidate
    
    os.makedirs(f"{work_dir}/{dir_name}", exist_ok=True)
    
    commands = []
    for file in fil:
        bandname = file.split('/')[-1].split('.')[0]
        # if p0 > 0.1, then use a -slow flag
        # Each command is submitted as a separate job
        commands.extend([
            f"prepfold -topo -dm {candidate_details['dm_opt']} -f {candidate_details['f0_opt']} -fd {candidate_details['f1_opt']} -o {dir_name}/{bandname}_{name} {file}",
            f"prepfold -topo -dm 0 -f {candidate_details['f0_opt']} -fd {candidate_details['f1_opt']} -o {dir_name}/{bandname}_{name}_0dm {file}",
            f"prepfold -topo -fine -dm {candidate_details['dm_opt']} -f {candidate_details['f0_opt']} -fd {candidate_details['f1_opt']} -o {dir_name}/{bandname}_{name}_fine {file}",
            f"prepfold -topo -nosearch -dm {candidate_details['dm_opt']} -f {candidate_details['f0_opt']} -fd {candidate_details['f1_opt']} -o {dir_name}/{bandname}_{name}_nosearch {file}",
            f"prepfold -topo -pfact 2 -dm {candidate_details['dm_opt']} -f {candidate_details['f0_opt']} -fd {candidate_details['f1_opt']} -o {dir_name}/{bandname}_{name}_pfact2 {file}",
            f"prepfold -topo -ffact 2 -dm {candidate_details['dm_opt']} -f {candidate_details['f0_opt']} -fd {candidate_details['f1_opt']} -o {dir_name}/{bandname}_{name}_ffact2 {file}"
        ])
        
    # num_cores = min(num_cores, len(commands)) # Use the minimum of the number of commands and the number of cores
    print(f"Number of cores: {num_cores}")
    # Use multiprocessing to submit multiple jobs at once
    pool = Pool(processes = num_cores)
    pool.map(submit_job, commands)
    pool.close()
    pool.join()
    
    end_time = time.time()
    print(f"Time taken: {end_time - start_time} seconds")
    

def main():
    parser = argparse.ArgumentParser(description="Candylicker creates bash script to check the legitimacy of a candidate. It's essentially licking to see if it is sweet.")
    parser.add_argument('-wd', type=str, help="Path to the work directory", default=os.getcwd())
    parser.add_argument('-csv', type=str, help="Path to the candidates.csv", default="")
    parser.add_argument('-cand', type=str, help="Path of the candidate png as per candyjar, the code uses this to search for pd and fd")
    parser.add_argument('-fil', type=str, help="Path to the fil files to fold, just give the path to the directory.")
    parser.add_argument('-cand_csv', type=str, help="Path to the csv file with candidates you want to fold obtained from candyjar", default="")
    parser.add_argument('-num_cores', type=int, help="Number of cores to use for multiprocessing", default=30)
    args = parser.parse_args()
    
    work_dir = args.wd
    csv_file = args.csv
    fil = glob.glob(args.fil + "/*.fil")
    cand = args.cand
    
    create_submit_file(work_dir, cand, csv_file, fil, args.num_cores)

if __name__ == "__main__":
    main()