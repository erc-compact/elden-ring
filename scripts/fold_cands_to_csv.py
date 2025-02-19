import argparse
import pandas as pd

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

def pulsarx_to_csv_format(filenames, pulsarx_cand_file, output_file):
    df = pd.read_csv(pulsarx_cand_file, skiprows=11, sep='\s+')
    df['fold_cands_filename'] = filenames

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
        "p0_new": float, "p1_new": float
    })
    
    df.to_csv(output_file, index=False)

def main():
    parser = argparse.ArgumentParser(description="Process pulsar data")
    # parser.add_argument("-p", "--process_name", choices=["pulsarx", "presto"], required=True, help="Name of the process to execute")
    parser.add_argument("-f", "--filenames", nargs='+', required=True, help="List of filenames to process")
    parser.add_argument("-c", "--pulsarx_cand_file", help="PulsarX candidate file")
    parser.add_argument("-o", "--output_file", help="Output file name", default="search_fold_cands.csv")

    args = parser.parse_args()
    
    pulsarx_to_csv_format(
            args.filenames, args.pulsarx_cand_file, args.output_file
        )
    
if __name__ == "__main__":
    main()