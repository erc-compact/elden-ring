import pandas as pd
import xml.etree.ElementTree as ET
from astropy import units as u
from astropy.coordinates import SkyCoord
import argparse
import os
import json


def read_candidate_files(files, chunk_id, verbose=False):
    """
    Reads XML candidate files and aggragates the candidates in a single pandas
        DataFrame.

    Args:
        files (list): List of file paths to read.
        verbose (bool): If True, prints additional debug information.

    Returns:
        df_candidates (DataFrame): DataFrame containing all candidates from the files.
        obs_meta_data (dict): Dictionary containing metadata.
    """

    # Number of files found
    if verbose:
        print(f"{len(files)} candidates files found.")

    # Initialize list to store all rows and counters
    all_rows = []
    file_index = 0
    corrupted_list = []
    obs_meta_data = {"tsamp": None,
                    "nsamples": None,
                    "obs_length": None,
                    "fft_size": None,           
                    'obs_length_over_c': None}

    # Loop over all files
    for file in files:
        file = file.replace(',','') 
        try:
            if verbose:
                print("Reading file: {}".format(file))
            # Parse XML file
            tree = ET.parse(file)
        except:
            # If file is corrupted, skip it and increment corrupted count
            if verbose:
                print("{} is corrupted, skipping file".format(file))
            corrupted_list.append(file)
            continue
        root = tree.getroot()
        header_params = root[1]
        search_params = root[2]
        segment_params = root[3]

        # Read the fil file from the XML file
        fil_file = root[2].find("infilename").text

        # Get candidates from the XML file
        candidates = root[7]
        
        # Create a row for each candidate and add it to the total list
        all_rows.extend(create_row(root, candidates, file, file_index, fil_file))
        

        # Grab needed meta data
        if file_index == 0:
            # Get metadata from the XML file
            tsamp = float(root[1].find("tsamp").text)
            total_nsamples = float(root[1].find("nsamples").text)
            segment_start_sample = int(segment_params.find('segment_start_sample').text)
            segment_nsamples = int(segment_params.find('segment_nsamples').text)
            xml_segment_pepoch = float(segment_params.find('segment_pepoch').text)
            tstart = float(root[1].find("tstart").text)
            fft_size = float(root.find('search_parameters/size').text)
            obs_length = tsamp * total_nsamples
            speed_of_light = 299792458.0
            obs_length_over_c = obs_length / speed_of_light
            source_name = root[1].find("source_name").text
            obs_meta_data = {"tsamp": tsamp,
                             "segment_start_sample": segment_start_sample,
                             "segment_nsamples": segment_nsamples,
                             "xml_segment_pepoch": xml_segment_pepoch,
                             "chunk_id": chunk_id,
                             "total_nsamples": total_nsamples,
                             "obs_length": obs_length,
                             "fft_size": fft_size,
                             "tstart": tstart,
                             "source_name": source_name,      
                             'obs_length_over_c': obs_length_over_c}
        file_index += 1
    
    if verbose:
        corrupted_count = len(corrupted_list)
        print("{} inputs were corrupted ({}%).".format(
        corrupted_count, 100.0 * corrupted_count / len(files)))
    obs_meta_data['xml_files'] = [x for x in files if x not in corrupted_list]
    obs_meta_data["corrupted_xml_files"] = corrupted_list

    if len(all_rows) == 0:
        # If all files were corrupted, return an empty DataFrame
        return pd.DataFrame(all_rows), obs_meta_data

    # Create a DataFrame from the list of rows
    df_candidates = pd.DataFrame(all_rows)
    

    # Cast columns to the appropriate data types
    df_candidates = df_candidates.astype({"snr": float, "dm": float, "period": float,
                                          "acc": float, "nassoc": int})

    if verbose:
        print(f"{len(df_candidates)} candidates read.")

    # Sort DataFrame by snr
    df_candidates.sort_values('snr', inplace=True, ascending=False)
    df_candidates.reset_index(inplace=True, drop=True)

    return df_candidates, obs_meta_data


def create_row(root, candidates, file, file_index, fil_file):
    # Read a candidate file and creates data rows

    src_raj = float(root[1].find("src_raj").text)
    src_dej = float(root[1].find("src_dej").text)
    src_rajd, src_dejd = convert_to_deg(src_raj, src_dej)
    rows = []

    # Enter attributes that should be ignored here
    ignored_entries = ['candidate']
    #ignored_entries = ['candidate', 'byte_offset', 'opt_period', 'folded_snr']
    for candidate in candidates:
        new_dict = {}
        for can_entry in candidate.iter():
            if not can_entry.tag in ignored_entries:
                new_dict[can_entry.tag] = can_entry.text
        cand_id = candidate.attrib.get("id")
        new_dict['cand_id_in_file'] = cand_id
        new_dict['src_raj'] = src_raj
        new_dict['src_rajd'] = src_rajd
        new_dict['src_dej'] = src_dej
        new_dict['src_dejd'] = src_dejd
        new_dict['file_index'] = file_index
        new_dict['file'] = file
        new_dict["fil_file"] = fil_file
        rows.append(new_dict)
    return rows

def convert_to_deg(ra, dec):
    # Convert hour angle strings to degrees
    ra_deg = int(ra/10000)
    ra_min = int(ra/100) - 100*ra_deg
    ra_sec = ra - 10000*ra_deg - 100*ra_min
       
    dec_deg = int(dec/10000)
    dec_min = int(dec/100) - 100*dec_deg
    dec_sec = dec - 10000*dec_deg - 100*dec_min

    ra_coord = "{}:{}:{:.2f}".format(ra_deg, abs(ra_min), abs(ra_sec))
    #print(ra_coord)
    dec_coord = "{}:{}:{:.2f}".format(dec_deg, abs(dec_min), abs(dec_sec))
    #print(dec_coord)
   
    #coord = SkyCoord(ra_string + ' ' + dec_string, unit=(u.hourangle, u.deg))
    coord = SkyCoord('%s %s'%(ra_coord,dec_coord), unit=(u.hourangle, u.deg))
    return coord.ra.deg, coord.dec.deg

def main():
    parser = argparse.ArgumentParser(description='Read candidate files and aggregate them in a single pandas DataFrame.')
    parser.add_argument('files', metavar='files', type=str, nargs='+',
                        help='List of XML candidate files to read.')
    parser.add_argument('--verbose', action='store_true',
                        help='Print additional debug information.')
    parser.add_argument('--outdir', type=str, default=os.getcwd(),
                        help='Output directory for the aggregated DataFrame')
    parser.add_argument('--outfile', type=str, default='candidates.csv',
                        help='Name of the output file')
    parser.add_argument('--metafile', type=str, default='metafile.meta', 
                        help='Name of the metadata file')
    parser.add_argument('--chunk_id', type=str, default='10')
    args = parser.parse_args()
    cands, meta = read_candidate_files(args.files, args.chunk_id, args.verbose)
    
    chunk_id = args.chunk_id

    # Save the DataFrame to a CSV file
    outname = os.path.join(args.outdir, "ck" + str(chunk_id) + args.outfile)
    cands.to_csv(outname, index=False)
    
    # Save the metadata to a text file
    meta_outname = os.path.join(args.outdir, str(chunk_id) + args.metafile)
    json.dump(meta, open(meta_outname, 'w'))

if __name__ == "__main__":
    main()