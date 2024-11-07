import pandas as pd
import argparse
import os
import fnmatch
import shutil
from glob import glob
import tarfile


def get_isot_from_mjd(mjd):
    return Time(mjd, format='mjd', scale='utc').isot

def get_utc_from_filterbank_name(filterbank_name):
    datepart = filterbank_name.split("_")[1]
    return datepart[:10] + datepart[10].replace("-", "T") + datepart[11:]

def parse_cands(input_filename):
    with open(input_filename, 'r') as file:
        lines = file.readlines()
    
    # Extract header information
    header_info = {}
    data_start_index = 0
    for i, line in enumerate(lines):
        if line.startswith('#id'):
            data_start_index = i+1
            break
        elif line.startswith('#'):
            key, value = line[1:].strip().split(maxsplit=1)
            header_info[key] = value

    # Extract column names and data
    columns = lines[data_start_index - 1].strip().split()
    data = []
    for line in lines[data_start_index:]:
        data.append(line.strip().split())
    
    # Create DataFrame for the data
    df = pd.DataFrame(data, columns=columns)

    # Add header information as columns
    for key, value in header_info.items():
        df[key] = value

    #rename columns
    new_cols = {"#id": "id", "dm_old": "dm_user", "dm_new": "dm_opt",
                "dm_err": "dm_opt_err", "f0_err": "f0_opt_err", "f1_err": "f1_opt_err",
                "acc_err": "acc_opt_err",
                "f0_old": "f0_user", "f0_new": "f0_opt", "f1_old": "f1_user",
                "f1_new": "f1_opt", "acc_old": "acc_user", "acc_new": "acc_opt",
                "S/N": "sn_fft", "S/N_new": "sn_fold",
                "Filename": "filterbank_path", "Telescope": "telescope",
                "Source_name": "source_name", "Beam": "beam_name",
                "Date": "mjd_start", "RA": "ra", "DEC": "dec", "GL": "gl",
                "GB": "gb", "MaxDM_YMW16": "maxdm_ymw16", "Pepoch": "pepoch"}
    df.rename(columns=new_cols, inplace=True)
    #calculate isot

    df['utc_start'] = df['filterbank_path'].apply(get_utc_from_filterbank_name)
    # df['utc_start'] = df['mjd_start'].apply(get_isot_from_mjd)

    #format columns
    df['id'] = df['id'].astype(int)

    # get the beam id
    df["beam_name"] = os.path.basename(input_filename).split('_')[1]
    df["beam_id"] = int(os.path.basename(input_filename).split('_')[1][4:])

    return df


def get_pics_vals(pics_file):
    pics_df = pd.read_csv(pics_file)
    pics_df['id'] = pics_df['filename'].apply(lambda x: int(os.path.basename(x).split('_')[-1].split('.')[0]))
    pics_df['beam_name'] = pics_df['filename'].apply(lambda x: os.path.basename(x).split('_')[-2])
    # rename
    new_cols = {"clfl2_trapum_Ter5": "pics_trapum_ter5",
                'PALFA_MeerKAT_L_SBAND_Best_Fscore': 'pics_palfa_meerkat_l_sband_best_fscore',
                'clfl2_PALFA': 'pics_palfa',
                'MeerKAT_L_SBAND_COMBINED_Best_Recall': 'pics_meerkat_l_sband_combined_best_recall'}
    pics_df.rename(columns=new_cols, inplace=True)
    
    return pics_df


def find_pngs(inpath):
    pngs = glob(os.path.join(inpath, '*.png'))
    
    vals = []
    for png in pngs:
        id = int(os.path.basename(png).split('_')[-1].split('.')[0])
        cfbf = os.path.basename(png).split('_')[-2]
        vals.append((id, cfbf, png))

    # Change abs path to relative path
    # pngs = [os.path.basename(png, inpath) for png in pngs]

    pngs_df = pd.DataFrame(vals, columns=['id', 'beam_name', 'png_path'])
    return pngs_df


def create_tarball(output_filename, meta_path, png_files, tarball_name):
    with tarfile.open(tarball_name, "w:gz", dereference=True) as tar:
        # Add output CSV file
        tar.add(output_filename, arcname=os.path.basename(output_filename))
        
        # Add meta file
        tar.add(meta_path, arcname='metafiles/' + os.path.basename(meta_path))

        # Add PNG files
        for png in png_files:
            tar.add(png, arcname='plots/' + os.path.basename(png))


def find_directories_with_prefix(path, prefix):
    directories = []
    
    # Walk through the directory tree
    for root, dirs, files in os.walk(path):
        # Filter directories that start with the given prefix
        for dir_name in fnmatch.filter(dirs, f'{prefix}*'):
            directories.append(os.path.join(root, dir_name))
    
    return directories

def process_dir(input_dir, output_dir, beam_name=None):
    # add DM to the output directory


    # Find all .ar and .png files and copy them to the output directory
    ar_files = glob(os.path.join(input_dir, '*.ar'))
    png_files = glob(os.path.join(input_dir, '*.png'))

    # for ar_file in ar_files:
    #     new_ar_file = os.path.join(output_dir, f"{DM}_" + os.path.basename(ar_file))
    #     shutil.copy(ar_file, new_ar_file)


    # Find cands and pics in the original directory
    input_candfiles = glob(os.path.join(input_dir, '*.cands'))
    input_pics_filenames = glob(os.path.join(input_dir, '*_scored.csv'))
    
    pics_dfs = []
    for pics_file in input_pics_filenames:
        df = get_pics_vals(pics_file)
        pics_dfs.append(df)
    pics_df = pd.concat(pics_dfs)
    pics_df = pics_df.drop(columns=['filename'])

    cands_dfs = []
    for candfile in input_candfiles:
        cands_dfs.append(parse_cands(candfile))
    cands_df = pd.concat(cands_dfs)
    
    # Merge the two dataframes  
    merged_df = cands_df.merge(pics_df, on=['id', 'beam_name'], how='inner')

    # calculate png path
    png_df = find_pngs(input_dir)
    merged_df = merged_df.merge(png_df, on=['id', 'beam_name'], how='inner')
    # print(merged_df.columns)
    # Get the manual values

    for png_file in png_files:
        new_png_file = os.path.join(output_dir, os.path.basename(png_file))
        shutil.copy(png_file, new_png_file)

    merged_df['pointing_id'] = args.pointing_id


    required_columns = ["pointing_id", "beam_id", "beam_name", "source_name",
                        "ra", "dec", "gl", "gb", "mjd_start", "utc_start",
                        "f0_user", "f0_opt", "f0_opt_err", "f1_user", "f1_opt",
                        "f1_opt_err", "acc_user", "acc_opt", "acc_opt_err",
                        "dm_user", "dm_opt", "dm_opt_err", "sn_fft", "sn_fold",
                        "pepoch", "maxdm_ymw16", "dist_ymw16", "pics_trapum_ter5",
                        "pics_palfa", "pics_palfa_meerkat_l_sband_best_fscore", "pics_meerkat_l_sband_combined_best_recall",
                        "png_path", "metafile_path",
                        "filterbank_path", "candidate_tarball_path"]

    # for name in required_columns:
    #     if name not in merged_df.columns:
    #         print(f"Column {name} not found in the dataframe")

    # Rename the png files
    try:
        merged_df.dropna(subset=['png_path'], inplace=True)
        merged_df['png_path'] = merged_df['png_path'].apply(lambda x: output_dir + '/' + os.path.basename(x))
    except:
        print(input_dir)

    # # order the columns
    # try:
    #     breakpoint()
    #     merged_df = merged_df[required_columns]
    # except:
    #     print(args.dir)

    return merged_df

def main(args):
    input_dir = args.dir
    output_dir = args.outdir

    # Create outdir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    total_df = process_dir(input_dir, output_dir)

    # Collect all files to be added to the tarball
    png_files = list(total_df['png_path'])
    meta_path = args.meta_path

    # Change png paths and meta paths to relative paths
    total_df.dropna(subset=['png_path'], inplace=True)
    total_df['png_path'] = total_df['png_path'].apply(lambda x: "plots/" + os.path.basename(x))
    total_df['metafile_path'] = "metafiles/" + os.path.basename(args.meta_path)

    DM_str = input_dir.split('/')[-1]
    beam_name = input_dir.split('/')[-2]
    
    tarball_name = os.path.join(output_dir, args.name)

    total_df['candidate_tarball_path'] = tarball_name

    required_columns = ["pointing_id", "beam_id", "beam_name", "source_name",
                    "ra", "dec", "gl", "gb", "mjd_start", "utc_start",
                    "f0_user", "f0_opt", "f0_opt_err", "f1_user", "f1_opt",
                    "f1_opt_err", "acc_user", "acc_opt", "acc_opt_err",
                    "dm_user", "dm_opt", "dm_opt_err", "sn_fft", "sn_fold",
                    "pepoch", "maxdm_ymw16", "dist_ymw16", "pics_trapum_ter5",
                    "pics_palfa", "pics_palfa_meerkat_l_sband_best_fscore", "pics_meerkat_l_sband_combined_best_recall",
                    "png_path", "metafile_path",
                    "filterbank_path", "candidate_tarball_path"]
    
    #order the columns
    total_df = total_df[required_columns]

    # save the final csv
    output_filename = os.path.join(output_dir, 'candidates.csv')
    total_df.to_csv(output_filename, index=False, header=True)

    # Create the tarball
    create_tarball(output_filename, meta_path, png_files, tarball_name)

arg_parser = argparse.ArgumentParser(
    description="A utility tool to obtain a csv file from presto candidates to view with CandyJar ")
arg_parser.add_argument("-d","--dir", help="File directory", required=True)
arg_parser.add_argument("-m", "--meta", dest="meta_path",
                        help="apsuse.meta file of that observation", required=True)
arg_parser.add_argument("-p", "--pointing_id", help="Pointing ID", required=True)
arg_parser.add_argument("-o", "--outdir", help="Output tarball",
                        default=os.getcwd()) 
arg_parser.add_argument("-n", "--name", help="Name of the output tarball",
                        default="candidates.tar.gz")
arg_parser.add_argument("-v", "--verbose", action="store_true",help="Verbose output")
arg_parser.add_argument("-g", '--global_beams', action='store_true',
                        default=False,
                        help='Create global beam candidates csv file')
args = arg_parser.parse_args()

if __name__ == "__main__":
    main(args)