#!/usr/bin/env python2.7

import os
import sys
import glob
import cPickle
import argparse
import pandas as pd
import logging
import subprocess
sys.path.append('/home/psr')
from ubc_AI.data import pfdreader
import errno
import theano


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def run_pics_parallel(filenames, pics_model, model_name):
    logging.debug("Loading model from %s", pics_model)
    classifier_model = cPickle.load(open(pics_model, 'rb'))
    try:
        logging.debug("Scoring files in parallel for model: %s", model_name)
        AI_scores = classifier_model.report_score([pfdreader(f) for f in filenames])
        df = pd.DataFrame({'filename': filenames, model_name: AI_scores})
        return df
    except Exception as e:
        logging.error("run_pics_parallel failed for model %s. Reason: %s", model_name, e)
        raise

def run_pics_sequential(filenames, pics_model, model_name):
    logging.debug("Loading model from %s", pics_model)
    classifier_model = cPickle.load(open(pics_model, 'rb'))
    AI_scores = []
    logging.debug("Scoring files sequentially for model: %s", model_name)
    for f in filenames:
        try:
            score = classifier_model.report_score([pfdreader(f)])
            AI_scores.append(score)
        except Exception as e:
            logging.error("Scoring failed for file %s. Reason: %s", f, e)
            AI_scores.append(0)
    df = pd.DataFrame({'filename': filenames, model_name: AI_scores})
    return df

def main(input_dir, pics_model_dir, fold_csv, fold_cands_col, generate_uuid, create_symlinks, sourcedir, publishdir, verbose):
    logging.info("Searching for input files (.pfd, .ar, .ar2) in %s", input_dir)
    filenames = (
        glob.glob(os.path.join(input_dir, '*.pfd')) +
        glob.glob(os.path.join(input_dir, '*.ar')) +
        glob.glob(os.path.join(input_dir, '*.ar2'))
    )
    logging.debug("Found %d candidate files.", len(filenames))

    # This DataFrame will hold the combined scores from all models
    master_df = pd.DataFrame({'filename': filenames})

    logging.info("Looking for PICS models in %s", pics_model_dir)
    models = glob.glob(os.path.join(pics_model_dir, '*.pkl'))
    logging.debug("Found models: %s", models)

    # Score all models and merge them horizontally into master_df
    for pics_model in models:
        model_name = os.path.splitext(os.path.basename(pics_model))[0]
        logging.info("Scoring with model: %s", model_name)
        try:
            df = run_pics_parallel(filenames, pics_model, model_name)
            logging.info("Parallel scoring successful for model: %s", model_name)
        except Exception:
            logging.info("Parallel scoring failed, attempting sequential scoring for model: %s", model_name)
            df = run_pics_sequential(filenames, pics_model, model_name)

        # If user requested UUID generation, create a new column for each model
        if generate_uuid:
            uuid_col_name = "cand_tracker_database_uuid_{}".format(model_name)
            df[uuid_col_name] = UUIDUtility.generate_uuid_list(len(df))
            #df[uuid_col_name] = [UUIDUtility.generate_uuid_string() for _ in range(len(df))]
        
        #Change filename to basename
        df['filename'] = df['filename'].apply(os.path.basename)
        master_df['filename'] = master_df['filename'].apply(os.path.basename)
        master_df = pd.merge(master_df, df, on='filename', how='left')
        logging.info("Merged model %s into master_df." % model_name)
        # If user requested symlinks, create them in publishdir
        if create_symlinks:
            logging.info("Creating symlinks to publishdir")
            try:
                model_symlink_dir = os.path.join(publishdir, model_name)
                mkdir_p(model_symlink_dir)
                shortlist = df.loc[df[model_name] >= 0.1]
                short_list_filenames = shortlist['filename'].values
                png_filenames = [os.path.splitext(f)[0] + '.png' for f in short_list_filenames]
                for f in png_filenames:
                    try:
                        source = os.path.join(sourcedir, f)
                        destination = os.path.join(model_symlink_dir, f)
                        #Skip if symlink already exists
                        if os.path.exists(destination):
                            logging.info("Symlink already exists for file %s", f)
                            continue
                        os.symlink(source, destination)
                    except Exception as e:
                        logging.error("Symlink creation failed for file %s. Reason: %s", f, e)
                        continue
            except Exception as e:
                logging.error("Symlink creation failed for model %s. Reason: %s", model_name, e)
                continue
   
    # If user provided fold_csv, do an *inner* merge with it, but invert the order:
    #   merged = pd.merge(merge_df, master_df, left_on=fold_cands_col, right_on='filename', how='inner')
    if fold_csv:
        logging.info("Merging with additional CSV: %s using column '%s'", fold_csv, fold_cands_col)
        try:
            merge_df = pd.read_csv(fold_csv)
            merged = pd.merge(
                merge_df, 
                master_df, 
                left_on=fold_cands_col, 
                right_on='filename', 
                how='inner'
            )
            master_df = merged
            logging.info("Merge successful. Final row count: %d", len(master_df))
        except Exception as e:
            logging.error("Merge failed: %s", e)
            raise

    # Decide on the output filename
    output_path = 'search_fold_pics_merged.csv' if fold_csv else 'pics_scores.csv'
    master_df.to_csv(output_path, index=False)
    logging.info("Scores have been written to: %s", os.path.abspath(output_path))

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run PICS scoring with optional merging and UUID generation.')
    parser.add_argument('-i', '--input_dir', required=True, help='Input directory of candidate files (.pfd, .ar, .ar2)')
    parser.add_argument('-m', '--model_dir', required=True, help='Directory containing PICS model files (*.pkl)')
    parser.add_argument('-f', '--fold_csv', default=None, 
                        help='Optional CSV to merge pics score with. If omitted, no merge is done.')
    parser.add_argument('--fold_cands_col', default='fold_cands_filename', 
                        help='Column name in the fold CSV to join on (default: fold_cands_filename)')
    parser.add_argument('-g', '--generate_uuid', action='store_true', 
                        help='Generate a unique UUID string for each candidate in each model.')
    parser.add_argument('-c', '--create_symlinks', action='store_true', help='Create symlinks to publishdir')
    parser.add_argument('-s', '--sourcedir', help='Source Directory for symlinks', default=os.getcwd())
    parser.add_argument('-p', '--publishdir', help='Directory to create symlinks in', default='publishdir')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose logging.')
   
    args = parser.parse_args()

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(levelname)s: %(message)s')

    main(args.input_dir, args.model_dir, args.fold_csv, args.fold_cands_col, args.generate_uuid, args.create_symlinks, args.sourcedir, args.publishdir, args.verbose)
