import argparse
import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET
from uncertainties import ufloat
from uncertainties.unumpy import uarray, nominal_values, std_devs

def a_to_pdot(P_s, acc_ms2):
    LIGHT_SPEED = 2.99792458e8  # Speed of Light in SI
    return P_s * acc_ms2 / LIGHT_SPEED

def calculate_spin(f=None, fdot=None, p=None, pdot=None):
    if f is not None and fdot is not None:
        p = 1 / f
        pdot = -fdot / (f ** 2)
    elif p is not None and pdot is not None:
        f = 1 / p
        fdot = -pdot * (p ** 2)
    else:
        raise ValueError("Either (f, fdot) or (p, pdot) must be provided")
    return f, fdot, p, pdot

def calculate_spin_with_error(f=None, fdot=None, p=None, pdot=None, f0_err=None, f1_err=None, p_err=None, pdot_err=None):
    """
    Compute spin parameters with optional Gaussian error propagation.
    
    Supports both scalars and NumPy arrays.

    Parameters:
    - f, fdot: Frequency and frequency derivative.
    - p, pdot: Period and period derivative.
    - f0_err, f1_err, p_err, pdot_err: Corresponding errors.

    Returns:
    - f, fdot, p, pdot (values)
    - f_err, fdot_err, p_err, pdot_err (errors, if applicable)
    """

    if f is not None and fdot is not None:
        f, fdot = np.asarray(f), np.asarray(fdot)  # Ensure NumPy arrays
        if f0_err is not None and f1_err is not None:
            f0_err, f1_err = np.asarray(f0_err), np.asarray(f1_err)
            f_u = uarray(f, f0_err)
            fdot_u = uarray(fdot, f1_err)

            p_u = 1 / f_u
            pdot_u = -fdot_u / (f_u ** 2)

            return nominal_values(p_u), nominal_values(pdot_u), std_devs(p_u), std_devs(pdot_u)

        else:
            p = 1 / f
            pdot = -fdot / (f ** 2)
            return p, pdot, None, None

    elif p is not None and pdot is not None:
        p, pdot = np.asarray(p), np.asarray(pdot)
        if p_err is not None and pdot_err is not None:
            p_err, pdot_err = np.asarray(p_err), np.asarray(pdot_err)
            p_u = uarray(p, p_err)
            pdot_u = uarray(pdot, pdot_err)

            f_u = 1 / p_u
            fdot_u = -pdot_u / (p_u ** 2)

            return nominal_values(f_u), nominal_values(fdot_u), std_devs(f_u), std_devs(fdot_u)

        else:
            f = 1 / p
            fdot = -pdot / (p ** 2)
            return f, fdot, None, None

    else:
        raise ValueError("Either (f, fdot) or (p, pdot) must be provided")
    

def pulsarx_to_csv_format(filenames, df, search_candidates_file, output_file, publish_dir):
    fold_data = df
    fold_data['fold_cands_filename'] = filenames
    fold_data['fold_cands_filepath'] = publish_dir
    search_data = pd.read_csv(search_candidates_file)
    search_data['index'] = search_data.index
    df = pd.merge(fold_data, search_data, left_on='#id', right_on='index', how='inner')
    del df['index']
    
    p, pdot, p_error, pdot_error = calculate_spin_with_error(f=df['f0_new'].values, fdot=df['f1_new'].values, f0_err=df['f0_err'].values, f1_err=df['f1_err'].values)

    df['p0_new'] = p
    df['p1_new'] = pdot
    df['p0_err'] = p_error
    df['p1_err'] = pdot_error
    
    f, fdot, p, pdot = calculate_spin(f=df['f0_new'].values, fdot=df['f1_new'].values)
    df['p0_new'] = p
    df['p1_new'] = pdot
    f, fdot, p, pdot = calculate_spin(f=df['f0_old'].values, fdot=df['f1_old'].values)
    df['p0_old'] = p
    df['p1_old'] = pdot
    df = df.astype({
        "#id": int, "fold_cands_filename": str, "f0_new": float, "f1_new": float,
        "dm_new": float, "S/N_new": float, "f0_old": float, "f1_old": float,
        "dm_old": float, "S/N": float, "p0_old": float, "p1_old": float,
        "p0_new": float, "p1_new": float, "p0_err" : float, "p1_err" : float
    })
    
    df.to_csv(output_file, index=False)

def main():
    parser = argparse.ArgumentParser(description="Process pulsar data")
    parser.add_argument("-f", "--filenames", nargs='+', required=True, help="List of filenames to process")
    parser.add_argument("-c", "--pulsarx_cand_files", nargs='+', help="PulsarX candidate files")
    parser.add_argument("-x", "--filtered_search_csv", required=True, help="Search CSV candidates selected from XML for folding")
    parser.add_argument("-o", "--output_file", help="Output file name", default="search_fold_cands.csv")
    parser.add_argument("-p", "--publish_dir", help="PNG directory", required=True)

    args = parser.parse_args()
        
    # create a master candidate file
    if args.pulsarx_cand_files:
        master_basename = f"{'.'.join(args.pulsarx_cand_files[0].split('/')[-1].split('.')[:2])}_master.cands"
        master_df = pd.DataFrame()
        current_id = 1
        for pulsarx_cand_file in args.pulsarx_cand_files:
            df = pd.read_csv(pulsarx_cand_file, skiprows=11, sep=r'\s+')
            df["#id"] = range(current_id, current_id + len(df))
            current_id += len(df)
            master_df = pd.concat([master_df, df], ignore_index=True)
        master_df.to_csv(f"{master_basename}", index=False)
    else:
        # error out
        print("No PulsarX candidate files provided. Please provide at least one.")
        
    pulsarx_to_csv_format(
            args.filenames, master_df, args.filtered_search_csv, args.output_file, args.publish_dir
        )
    
if __name__ == "__main__":
    main()