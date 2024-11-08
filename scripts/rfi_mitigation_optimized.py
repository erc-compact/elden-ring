import numpy as np
import numpy as np
import os
from itertools import islice

import struct
from astropy.io import fits
import matplotlib.pyplot as plt
from collections import Counter
from scipy.stats import kurtosis
import os
import argparse
from concurrent.futures import ProcessPoolExecutor
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
from astropy.time import Time
import astropy.units as u
from astropy.time import Time, TimeDelta
from sigpyproc.readers import FilReader

# Filterbank header parsing functions
def _read_int(fh):
    return struct.unpack('i', fh.read(4))[0]

def _read_double(fh):
    return struct.unpack('d', fh.read(8))[0]

def _read_string(fh):
    length = _read_int(fh)
    return fh.read(length)

def read_data_file(file_path, start_time, end_time, file_type='fits'):
    """
    Read data from either FITS or filterbank file
    """
    if file_type == 'fits':
        return read_fits_file(file_path, start_time, end_time)
    else:
        return read_filterbank_file(file_path, start_time, end_time)
        
# Read FITS File
def read_fits_file(fits_file, start_time, end_time):
    with fits.open(fits_file) as hdul:
        data = hdul['SUBINT'].data['DATA']
        header = hdul['SUBINT'].header
        nchan = header.get('NCHAN')

        nsubint = len(data)
        time_resolution_s = header['TBIN']
        time_resolution_us = time_resolution_s * 1e6
        time_per_subint = time_resolution_us * data.shape[1]

        #BAND3
        # min_freq_mhz = 2976
        # max_freq_mhz = 4101
        
        hdr = hdul[0].header
        OBSFREQ = hdr['OBSFREQ']
        OBSBW = abs(hdr['OBSBW'])
        min_freq_mhz = OBSFREQ - OBSBW/2
        max_freq_mhz = OBSFREQ + OBSBW/2        

        # min_freq_mhz = 1290
        # max_freq_mhz = 1940

        # Convert start_time and end_time from seconds to indices in the data array
        start_index = int(start_time * 1e6 / time_per_subint)
        end_index = int(end_time * 1e6 / time_per_subint)

        start_index = max(0, start_index)
        end_index = min(nsubint, end_index)

        data = data[start_index:end_index, :, 0]
        dynamic_spectrum = data.reshape((data.shape[0] * data.shape[1], nchan))

        time_bins_per_block = int(1 / time_resolution_s)
        processed_data_size_bytes = dynamic_spectrum.nbytes

        print(f"Processed data size: {processed_data_size_bytes} bytes")
        return dynamic_spectrum, time_bins_per_block, int(1e6 / time_resolution_us), nchan, min_freq_mhz, max_freq_mhz, time_resolution_us


def read_filterbank_header(filename):
    """Read filterbank header and return header info as a dictionary."""
    header = {}
    with open(filename, 'rb') as fh:
        # Check if file starts with 'HEADER_START'
        keyword = _read_string(fh)
        if keyword != 'HEADER_START':
            raise ValueError("Not a valid filterbank file")

        while True:
            keyword = _read_string(fh)
            if keyword == 'HEADER_END':
                break

            if keyword in ['source_name', 'rawdatafile']:
                header[keyword] = _read_string(fh)
            elif keyword in ['az_start', 'za_start', 'src_raj', 'src_dej', 'tstart', 'tsamp',
                           'fch1', 'foff', 'nchans', 'nifs', 'nbits', 'nbeams',
                           'ibeam', 'machine_id', 'telescope_id', 'data_type']:
                header[keyword] = _read_double(fh)

        # Calculate header size
        header['header_size'] = fh.tell()

    return header

def get_file_header(file_path):
    """
    Get header information from either FITS or filterbank file
    """
    if file_path.endswith('.fil'):
        fil = FilReader(file_path)
        # Convert telescope coordinates to RA/DEC if available
        try:
            ra = fil.header.ra_deg if hasattr(fil.header, 'ra_deg') else fil.header.src_ra
            dec = fil.header.dec_deg if hasattr(fil.header, 'dec_deg') else fil.header.src_dec
        except AttributeError:
            try:
                ra = fil.header.tstart_RA if hasattr(fil.header, 'tstart_RA') else 0.0
                dec = fil.header.tstart_DEC if hasattr(fil.header, 'tstart_DEC') else 0.0
            except AttributeError:
                ra = 0.0
                dec = 0.0
                print("Warning: Could not find RA/DEC information in filterbank header")

        header = {
            'STT_IMJD': int(fil.header.tstart),  # Integer MJD
            'STT_SMJD': int((fil.header.tstart % 1) * 86400),  # Seconds portion
            'STT_OFFS': ((fil.header.tstart % 1) * 86400) % 1,  # Fractional seconds
            'RA': ra,
            'DEC': dec,
            'TBIN': fil.header.tsamp,  # Sampling time
            'NAXIS2': fil.header.nsamples  # Number of samples
        }
        return header
    else:
        with fits.open(file_path) as hdul:
            primary_header = hdul[0].header
            subint_header = hdul['SUBINT'].header
            header = {
                'STT_IMJD': primary_header['STT_IMJD'],
                'STT_SMJD': primary_header['STT_SMJD'],
                'STT_OFFS': primary_header['STT_OFFS'],
                'RA': primary_header['RA'],
                'DEC': primary_header['DEC'],
                'TBIN': subint_header['TBIN'],
                'NAXIS2': subint_header['NAXIS2']
            }
            return header

def read_filterbank_file(filename, start_time, end_time):
    """
    Read filterbank data for a specific time range using sigpyproc.
    
    Parameters:
    filename (str): Path to filterbank file
    start_time (float): Start time in seconds
    end_time (float): End time in seconds
    
    Returns:
    tuple: (dynamic_spectrum, time_bins_per_block, sampling_rate, nchan, min_freq_mhz, max_freq_mhz, time_resolution_us)
    """
    # Initialize FilReader
    fil = FilReader(filename)
    
    # Get header information
    tsamp = fil.header.tsamp  # Sampling time in seconds
    nchan = fil.header.nchans  # Number of channels
    fch1 = fil.header.fch1    # First channel frequency
    foff = fil.header.foff    # Channel bandwidth
    
    # Calculate frequency range
    frequencies = np.arange(nchan) * foff + fch1  # MHz
    min_freq_mhz = np.min(frequencies)
    max_freq_mhz = np.max(frequencies)
    
    # Calculate start and end samples
    start_sample = int(start_time / tsamp)
    end_sample = int(end_time / tsamp)
    n_samples = end_sample - start_sample
    
    # Read data in chunks to manage memory
    chunk_size = 10000  # Number of time samples per chunk
    n_chunks = (n_samples + chunk_size - 1) // chunk_size
    
    # Pre-allocate array for full data - note the order of dimensions
    dynamic_spectrum = np.empty((n_samples, nchan), dtype=np.float32)
    
    # Read data in chunks
    for i in range(n_chunks):
        chunk_start = start_sample + i * chunk_size
        chunk_end = min(chunk_start + chunk_size, end_sample)
        chunk_size_actual = chunk_end - chunk_start
        
        # Read chunk using read_block method
        block = fil.read_block(chunk_start, chunk_size_actual)
        # Convert FilterbankBlock to numpy array and transpose if necessary
        chunk_data = block.data
        if chunk_data.shape[0] == nchan:  # If dimensions are swapped
            chunk_data = chunk_data.T
        
        # Store in the pre-allocated array
        start_idx = i * chunk_size
        end_idx = start_idx + chunk_size_actual
        dynamic_spectrum[start_idx:end_idx] = chunk_data
    
    # Calculate return values
    time_resolution_us = tsamp * 1e6
    time_bins_per_block = int(1 / tsamp)
    sampling_rate = int(1 / tsamp)
    
    # Print memory usage info
    processed_data_size_bytes = dynamic_spectrum.nbytes
    print(f"Processed data size: {processed_data_size_bytes} bytes")
    
    return (
        dynamic_spectrum,
        time_bins_per_block,
        sampling_rate,
        nchan,
        min_freq_mhz,
        max_freq_mhz,
        time_resolution_us
    )
  
def calculate_az_el(fits_header, observation_time):
    ra_str = fits_header['RA']  
    dec_str = fits_header['DEC']  

    # Effelsberg Telescope Coordinates
    latitude = 50.5247  # Latitude in degrees
    longitude = 6.8836  # Longitude in degrees
    elevation = 369  # Altitude in meters
    location = EarthLocation(lat=latitude * u.deg, lon=longitude * u.deg, height=elevation * u.m)

    # Create a SkyCoord object for the celestial coordinates, with explicit units
    sky_coord = SkyCoord(ra=ra_str, dec=dec_str, unit=(u.hourangle, u.deg), frame='icrs')

    # Define the AltAz frame for the observer's location at the given time
    altaz_frame = AltAz(obstime=observation_time, location=location)

    # Convert the sky coordinates to AltAz
    altaz = sky_coord.transform_to(altaz_frame)

    # Extract the azimuth and elevation
    azimuth = altaz.az.degree
    elevation = altaz.alt.degree

    return azimuth, elevation

# Calculate Statistics
def calculate_statistics(dynamic_spectrum, time_bins_per_block):
    n_blocks = dynamic_spectrum.shape[0] // time_bins_per_block

    reshaped_spectrum = dynamic_spectrum[:n_blocks * time_bins_per_block].reshape((n_blocks, time_bins_per_block, -1))

    mean = np.mean(reshaped_spectrum, axis=1)
    std_dev = np.std(reshaped_spectrum, axis=1)
    sk = kurtosis(reshaped_spectrum, axis=1, fisher=True, bias=False)
    time_bins = time_bins_per_block
    print(time_bins)
    sk_variance = (24 * time_bins * (time_bins - 1) ** 2) / ((time_bins - 3) * (time_bins - 2) * (time_bins + 3) * (time_bins + 5))
    sk_variance = np.full_like(sk, sk_variance)

    skewness = np.mean((reshaped_spectrum - mean[:, np.newaxis, :]) ** 3, axis=1) / std_dev ** 3

    percentiles = np.percentile(reshaped_spectrum, [25, 50, 75, 100], axis=1)
    q1, q2, q3, q4 = percentiles
    print(n_blocks)

    return reshaped_spectrum, mean, std_dev, sk, sk_variance, skewness, q1, q2, q3, q4

# Save Statistics to File
def save_statistics(statistics_file, mean, std_dev, sk, sk_variance, skewness, q1, q2, q3, q4, frequencies):
    n_blocks, n_channels = mean.shape
    os.makedirs(os.path.dirname(statistics_file), exist_ok=True)

    header = f"{'Block':>8} {'Channel':>8} {'Frequency':>10} {'Mean':>12} {'Std Dev':>12} {'Sk':>12} {'Sk Variance':>12} {'Skewness':>12} {'Q1':>12} {'Q2':>12} {'Q3':>12} {'Q4':>12}\n"
    stats = np.column_stack([
        np.repeat(np.arange(n_blocks), n_channels),
        np.tile(np.arange(n_channels), n_blocks),
        np.tile(frequencies, n_blocks),
        mean.ravel(),
        std_dev.ravel(),
        sk.ravel(),
        sk_variance.ravel(),
        skewness.ravel(),
        q1.ravel(),
        q2.ravel(),
        q3.ravel(),
        q4.ravel()
    ])

    fmt = "%8d %8d %10.2f " + "%12.4f " * 9
    np.savetxt(statistics_file, stats, fmt=fmt, header=header)

# Identify Outliers using IQR
def identify_outliers_iqr(sk, q1, q3):
    iqr = q3 - q1
    lower_bound = q1 - 7 * iqr
    upper_bound = q3 + 7 * iqr
    outliers_iqr = (sk < lower_bound) | (sk > upper_bound)
    return outliers_iqr

# Identify Outliers using Spectral Kurtosis (SK)
def identify_outliers_sk(sk, sk_variance):
    sk_threshold = 1 + 2.6 * sk_variance
    outliers_sk = np.abs(sk) > sk_threshold
    return outliers_sk

# Identify Frequent Outliers and Expand with Neighbors
def identify_frequent_outliers(outliers_sk, frequencies, threshold=0.7, neighbor_threshold=0.5):
    occurrence_count = np.sum(outliers_sk, axis=0)
    frequent_outliers = np.where(occurrence_count > (threshold * outliers_sk.shape[0]))[0]

    # Expand the frequent outliers by checking neighboring channels
    expanded_frequent_outliers = np.copy(frequent_outliers)

    left_neighbors = frequent_outliers - 1
    left_neighbors = left_neighbors[(left_neighbors >= 0)]
    left_neighbors = left_neighbors[(occurrence_count[left_neighbors] > (neighbor_threshold * outliers_sk.shape[0]))]

    right_neighbors = frequent_outliers + 1
    right_neighbors = right_neighbors[(right_neighbors < len(frequencies))]
    right_neighbors = right_neighbors[(occurrence_count[right_neighbors] > (neighbor_threshold * outliers_sk.shape[0]))]

    expanded_frequent_outliers = np.unique(np.concatenate((expanded_frequent_outliers, left_neighbors, right_neighbors)))

    # Convert the set of outliers to contiguous frequency ranges
    expanded_frequent_outliers.sort()
    diff = np.diff(expanded_frequent_outliers)
    split_indices = np.where(diff > 1)[0] + 1
    range_start_indices = np.insert(split_indices, 0, 0)
    range_end_indices = np.append(split_indices, len(expanded_frequent_outliers))

    frequency_ranges_in_mhz = [(frequencies[expanded_frequent_outliers[start]], frequencies[expanded_frequent_outliers[end - 1]])
                               for start, end in zip(range_start_indices, range_end_indices)]

    return frequency_ranges_in_mhz, expanded_frequent_outliers

# Function to Save Outliers to File
def save_outliers_to_file(outliers, filename, frequencies=None, is_freq=False):

    os.makedirs(os.path.dirname(filename), exist_ok=True)

    if is_freq:
        with open(filename, 'w') as f:
            f.write("Frequency Ranges (MHz):\n")
            for start, end in outliers:
                f.write(f"{start:.2f} - {end:.2f} MHz\n")
    else:
        block_indices, channel_indices = np.nonzero(outliers)
        outlier_frequencies = frequencies[channel_indices]
        output_data = np.column_stack((block_indices, outlier_frequencies))
        np.savetxt(filename, output_data, fmt="%d %.2f MHz", header="Block Frequency(MHz)", comments="")

# Merge Frequency Ranges
def merge_frequency_ranges(frequency_ranges):
    if not frequency_ranges:
        return []

    frequency_ranges = sorted(frequency_ranges, key=lambda x: x[0])
    merged_ranges = [frequency_ranges[0]]

    for current_start, current_end in frequency_ranges[1:]:
        last_start, last_end = merged_ranges[-1]
        if current_start <= last_end + 1:
            merged_ranges[-1] = (last_start, max(last_end, current_end))
        else:
            merged_ranges.append((current_start, current_end))

    return merged_ranges

# Plot Fraction of RFI Channels vs Time Blocks
def plot_fraction_rfi_vs_time(outliers_sk, output_folder):
    rfi_count_per_block = np.sum(outliers_sk, axis=1)
    time_blocks = np.arange(len(rfi_count_per_block))

    plt.figure(figsize=(10, 6))
    plt.scatter(time_blocks, rfi_count_per_block, color='blue', alpha=0.6)
    plt.xlabel('Time Block Index')
    plt.ylabel('Count of RFI Channels')
    plt.title('Count of RFI Channels vs Time Blocks')
    plt.grid(True)

    output_file = os.path.join(output_folder, 'count_rfi_vs_time.png')
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()

# Plot Dynamic Spectrum Block
def plot_dynamic_spectrum_block(
    dynamic_spectrum_block, block_idx, outliers_sk, output_folder,
    sample_time_us, min_freq_mhz, max_freq_mhz, nchan, sk_block,
    frequencies, q3_block, q1_block, frequent_outliers_sk
):
    time_bins, nchan = dynamic_spectrum_block.shape[1:]
    extent = [0, time_bins * sample_time_us / 1e6, min_freq_mhz, max_freq_mhz]

    # Create a figure with 3 subplots: 2 for dynamic spectra, 1 for SK vs. frequency
    fig, axs = plt.subplots(
        nrows=1, ncols=3, figsize=(18, 8),
        gridspec_kw={'width_ratios': [1, 1, 0.25]}, sharey=True
    )

    # Plot for SK-based outliers (dynamic spectrum 1)
    im1 = axs[0].imshow(
        dynamic_spectrum_block[block_idx].T, aspect='auto', cmap='viridis', origin='lower', extent=extent,
        vmin=np.nanquantile(dynamic_spectrum_block[block_idx], 0.02),
        vmax=np.nanquantile(dynamic_spectrum_block[block_idx], 0.98)
    )
    axs[0].set_title(f'Dynamic Spectrum - Block {block_idx} (SK Outliers)')
    axs[0].set_xlabel('Time (s)')
    axs[0].set_ylabel('Frequency (MHz)')
    sk_outliers_frequencies = frequencies[outliers_sk[block_idx]]
    axs[0].hlines(sk_outliers_frequencies, *extent[:2], colors='red', linestyles='--', linewidth=0.6)

    # Plot for combined outliers with frequent SK outliers and their neighbors highlighted (dynamic spectrum 2)
    axs[1].imshow(
        dynamic_spectrum_block[block_idx].T, aspect='auto', cmap='viridis', origin='lower', extent=extent,
        vmin=np.nanquantile(dynamic_spectrum_block[block_idx], 0.02),
        vmax=np.nanquantile(dynamic_spectrum_block[block_idx], 0.98)
    )
    axs[1].set_title(f'Dynamic Spectrum - Block {block_idx} (Frequent Outliers + Neighbors)')
    axs[1].set_xlabel('Time (s)')
    for start_freq, end_freq in frequent_outliers_sk:
        axs[1].axhspan(start_freq, end_freq, color='orange', alpha=0.5)

    # Plot SK vs Frequency in the last subplot (narrower)
    axs[2].plot(sk_block[block_idx], frequencies, 'b-', linewidth=1)
    axs[2].set_xlabel('SK (Kurtosis)', color='blue')
    axs[2].tick_params(axis='x', labelcolor='orange')
    axs[2].set_yticks([])  # Hide y-axis ticks
    
    left_limit = sk_block[block_idx].min()
    right_limit = sk_block[block_idx].max()
    
    # Check if the limits are NaN or inf and set them to None if they are
    if np.isnan(left_limit) or np.isinf(left_limit):
        left_limit = None
    if np.isnan(right_limit) or np.isinf(right_limit):
        right_limit = None
    
    # Set the x-axis limits
    axs[2].set_xlim(left=left_limit, right=right_limit)
    # axs[2].set_xlim(left=sk_block[block_idx].min(), right=sk_block[block_idx].max())

    # Plot the colorbar only once for all subplots
    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    fig.colorbar(im1, cax=cbar_ax, label='Intensity')

    fig.tight_layout(rect=[0, 0, 0.9, 1])
    os.makedirs(output_folder, exist_ok=True)
    output_file = os.path.join(output_folder, f"dynamic_spectrum_block_{block_idx}.png")
    plt.savefig(output_file)
    plt.close()

# Plot Fraction of RFI Channels vs Time Blocks
def plot_fraction_rfi_vs_time(outliers_sk, output_folder):
    rfi_count_per_block = np.sum(outliers_sk, axis=1)
    time_blocks = np.arange(len(rfi_count_per_block))

    plt.figure(figsize=(10, 6))
    plt.scatter(time_blocks, rfi_count_per_block, color='blue', alpha=0.6)
    plt.xlabel('Time Block Index')
    plt.ylabel('Count of RFI Channels')
    plt.title('Count of RFI Channels vs Time Blocks')
    plt.grid(True)

    output_file = os.path.join(output_folder, 'count_rfi_vs_time.png')
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()

def process_interval(interval_idx, file_path, output_folder, target_resolution_ms, interval_duration, total_intervals, azimuth_start, elevation_start):
    start_time = interval_idx * interval_duration
    end_time = (interval_idx + 1) * interval_duration

    print(f"Running analysis for interval {interval_idx + 1}/{total_intervals}: Start time = {start_time}, End time = {end_time}")

    # Determine file type and read data accordingly
    file_type = 'filterbank' if file_path.endswith('.fil') else 'fits'
    dynamic_spectrum, time_bins_per_block, sampling_rate, nchan, min_freq_mhz, max_freq_mhz, sample_time_us = read_data_file(
        file_path, start_time, end_time, file_type
    )

    # Rest of the processing remains the same
    bins_per_ms = int(1e3 / sample_time_us)
    bins_per_target = int(target_resolution_ms * bins_per_ms)

    n_bins, n_channels = dynamic_spectrum.shape
    n_blocks = n_bins // bins_per_target
    reshaped_spectrum = dynamic_spectrum[:n_blocks * bins_per_target, :].reshape((n_blocks, bins_per_target, n_channels))

    dynamic_spectrum = reshaped_spectrum.mean(axis=1)

    new_time_resolution_us = target_resolution_ms * 1e3
    time_bins_per_block = int(1e6 / new_time_resolution_us)

    reshaped_spectrum, mean, std_dev, sk, sk_variance, skewness, q1, q2, q3, q4 = calculate_statistics(dynamic_spectrum, time_bins_per_block)

    outliers_iqr = identify_outliers_iqr(sk, q1, q3)
    outliers_sk = identify_outliers_sk(sk, sk_variance)

    frequencies = np.linspace(min_freq_mhz, max_freq_mhz, n_channels)
    frequent_outliers_sk, expanded_frequent_outliers_indices = identify_frequent_outliers(outliers_sk, frequencies)

    interval_output_folder = os.path.join(output_folder, f'interval_{interval_idx + 1}')
    os.makedirs(interval_output_folder, exist_ok=True)

    # Save results and create plots
    STATISTICS_FILE = os.path.join(interval_output_folder, 'statistics.txt')
    save_statistics(STATISTICS_FILE, mean, std_dev, sk, sk_variance, skewness, q1, q2, q3, q4, frequencies)

    OUTLIERS_SK_FILE = os.path.join(interval_output_folder, 'outliers_sk.txt')
    save_outliers_to_file(outliers_sk, OUTLIERS_SK_FILE, frequencies=frequencies)

    FREQUENT_OUTLIERS_FILE = os.path.join(interval_output_folder, 'frequent_outliers.txt')
    save_outliers_to_file(frequent_outliers_sk, FREQUENT_OUTLIERS_FILE, is_freq=True)

    reshaped_blocks = reshaped_spectrum.reshape((-1, time_bins_per_block, n_channels))
    for block_idx in range(reshaped_blocks.shape[0]):
        plot_dynamic_spectrum_block(
            reshaped_blocks, block_idx, outliers_sk, interval_output_folder,
            target_resolution_ms * 1e3, min_freq_mhz, max_freq_mhz, n_channels, 
            sk, frequencies, q3, q1, frequent_outliers_sk
        )

    plot_combined_sk_heatmap_and_histogram_interval(
        sk, sk_variance, frequencies, outliers_sk, frequent_outliers_sk,
        interval_output_folder, azimuth_start, elevation_start
    )

    return sk, sk_variance, outliers_sk, frequent_outliers_sk, frequencies

def plot_combined_sk_heatmap_and_histogram(
    sk, sk_variance, frequencies, outliers_sk, frequent_outliers_sk, output_folder,
    azimuth_start, elevation_start, azimuth_end, elevation_end, fits_header, fits_file
):
    num_blocks, num_channels = sk.shape  # Define num_channels here

    # Flatten outliers_sk to get corresponding frequencies for histogram
    outliers_sk_flat = np.flatnonzero(outliers_sk)
    combined_outliers_sk = frequencies[outliers_sk_flat % num_channels]

    # Create subplots for the text panel, histogram, and heatmap, with unequal height and shared x-axis
    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(15, 12), gridspec_kw={'height_ratios': [1, 1, 3]}, sharex=True)

    # Set font to a monospace font similar to the one in the image
    font_properties = {'family': 'monospace', 'size': 12}

    # Split the observation info into two columns
    left_info = (
        f"Telescope: Effelsberg\n"
        f"Source name: {os.path.basename(fits_file)}\n"
        f"RA: {fits_header.get('RA' )}\n"
        f"Dec: {fits_header.get('DEC')}\n"
        f"Duration of blocks in ms: 100 \n"


    )

    right_info = (
        f"Azimuth Start (deg): {azimuth_start:.2f}\n"
        f"Elevation Start (deg): {elevation_start:.2f}\n"
        f"Azimuth End (deg): {azimuth_end:.2f}\n"
        f"Elevation End (deg): {elevation_end:.2f}\n" 
    )

    # Create two text boxes, one for each column of information
    axs[0].text(0.05, 0.5, left_info, verticalalignment='center', horizontalalignment='left', transform=axs[0].transAxes, fontdict=font_properties)
    axs[0].text(0.55, 0.5, right_info, verticalalignment='center', horizontalalignment='left', transform=axs[0].transAxes, fontdict=font_properties)
    axs[0].axis('off')  # Turn off the axis for the text panel

    # Plot Histogram of bad frequencies (outliers)
    frequency_count = Counter(combined_outliers_sk)
    frequencies_list = np.array(sorted(frequency_count.keys()))
    counts_list = np.array([frequency_count[freq] for freq in frequencies_list])
    axs[1].bar(frequencies_list, counts_list, color='black', alpha=0.7, label='Outliers')

    # Add bars for frequent outliers + neighbors
    frequent_mask = np.zeros_like(frequencies_list, dtype=bool)
    for start, end in frequent_outliers_sk:
        mask = (frequencies_list >= start) & (frequencies_list <= end)
        frequent_mask |= mask
    axs[1].bar(frequencies_list[frequent_mask], counts_list[frequent_mask], color='orange', alpha=0.7, label='Outliers with 70% occurrence + their neighbors')

    axs[1].set_ylabel('Blocks')
    axs[1].set_title('Histogram of Bad Frequencies')
    axs[1].legend()

    # Plot combined SK heatmap
    cmap = plt.get_cmap('cividis')
    cmap.set_over('orange')  # Highlight values above the threshold

    # Define extent for the heatmap
    extent = [frequencies[0], frequencies[-1], 0, num_blocks]
    cax = axs[2].imshow(sk, aspect='auto', cmap=cmap, origin='lower', extent=extent, interpolation='none',
                        vmin=np.nanquantile(sk, 0.02), vmax=np.nanquantile(sk, 0.98))

    # Identify time blocks with more than 30% of the total frequency channels flagged as outliers
    percent_outliers_per_block = np.sum(outliers_sk, axis=1) / num_channels * 100
    blocks_with_high_outliers = np.where(percent_outliers_per_block > 30)[0]

    # Define block_data as a list of tuples containing block numbers and their percentage of bad channels
    block_data = [(block, percent_outliers_per_block[block]) for block in blocks_with_high_outliers]

    # Add segments for each identified time block based on their percentage
    for block, percent in block_data:
        color = 'blue' if 30 <= percent <= 50 else 'red'
        axs[2].plot([frequencies[0] - 5, frequencies[0]], [block, block], color=color, lw=2)  # Left segment
        axs[2].plot([frequencies[-1], frequencies[-1] + 5], [block, block], color=color, lw=2)  # Right segment
        #axs[2].text(frequencies[-1] + 10, block, f'{block}',
                   #color=color, fontsize=6, verticalalignment='center')

    axs[2].set_title('Combined Spectral Kurtosis Heatmap')
    axs[2].set_ylabel('Time Blocks')
    axs[2].set_xlabel('Frequency (MHz)')

    # Add the colorbar showing actual SK values
    cbar = fig.colorbar(cax, ax=axs[2], orientation='horizontal', pad=0.1, aspect=40)

    # Add a legend for the color-coded segments
    axs[2].plot([], [], color='blue', lw=2, label='30% - 50% bad channels')
    axs[2].plot([], [], color='red', lw=2, label='>50% bad channels')
    axs[2].legend(loc='upper right', fontsize=8)

    # Save combined plots
    output_file = os.path.join(output_folder, 'combined_sk_heatmap_and_histogram.png')
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()

    # Save block data to a text file
    block_data_file = os.path.join(output_folder, 'block_bad_channel_percentages.txt')
    with open(block_data_file, 'w') as f:
        f.write('Block Number | Percentage of Bad Channels (%)\n')
        for block, percent in block_data:
            f.write(f'{block:8d} | {percent:6.2f}\n')

def plot_combined_sk_heatmap_and_histogram_interval(
    sk, sk_variance, frequencies, outliers_sk, frequent_outliers_sk, output_folder,
    azimuth, elevation
):
    num_blocks, num_channels = sk.shape

    # Flatten outliers_sk to get corresponding frequencies for histogram
    outliers_sk_flat = np.flatnonzero(outliers_sk)
    combined_outliers_sk = frequencies[outliers_sk_flat % num_channels]

    # Create subplots for the histogram and heatmap, with unequal height and shared x-axis
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(15, 10), gridspec_kw={'height_ratios': [1, 3]}, sharex=True)

    # Plot Histogram of bad frequencies (outliers)
    frequency_count = Counter(combined_outliers_sk)
    frequencies_list = np.array(sorted(frequency_count.keys()))
    counts_list = np.array([frequency_count[freq] for freq in frequencies_list])
    axs[0].bar(frequencies_list, counts_list, color='black', alpha=0.7, label='Outliers')

    # Add bars for frequent outliers + neighbors
    frequent_mask = np.zeros_like(frequencies_list, dtype=bool)
    for start, end in frequent_outliers_sk:
        mask = (frequencies_list >= start) & (frequencies_list <= end)
        frequent_mask |= mask
    axs[0].bar(frequencies_list[frequent_mask], counts_list[frequent_mask], color='orange', alpha=0.7, label='Outliers with 70% occurrence + their neighbors')

    axs[0].set_ylabel('Blocks')
    axs[0].set_title('Histogram of Bad Frequencies for Interval')
    axs[0].legend()

    # Plot combined SK heatmap
    cmap = plt.get_cmap('cividis')
    cmap.set_over('orange')  # Highlight values above the threshold

    # Define extent for the heatmap
    extent = [frequencies[0], frequencies[-1], 0, num_blocks]
    cax = axs[1].imshow(sk, aspect='auto', cmap=cmap, origin='lower', extent=extent, interpolation='none',
                        vmin=np.nanquantile(sk, 0.02), vmax=np.nanquantile(sk, 0.98))

    # Identify time blocks with more than 30% of the total frequency channels flagged as outliers
    percent_outliers_per_block = np.sum(outliers_sk, axis=1) / num_channels * 100
    blocks_with_high_outliers = np.where(percent_outliers_per_block > 30)[0]

    # Define block_data as a list of tuples containing block numbers and their percentage of bad channels
    block_data = [(block, percent_outliers_per_block[block]) for block in blocks_with_high_outliers]

    # Add segments for each identified time block based on their percentage
    for block, percent in block_data:
        color = 'blue' if 30 <= percent <= 50 else 'red'
        axs[1].plot([frequencies[0] - 5, frequencies[0]], [block, block], color=color, lw=2)  # Left segment
        axs[1].plot([frequencies[-1], frequencies[-1] + 5], [block, block], color=color, lw=2)  # Right segment
        axs[1].text(frequencies[-1] + 10, block, f'{block}\n{percent:.1f}%',
                    color=color, fontsize=8, verticalalignment='center')

    # Set labels and titles
    axs[1].set_title('Spectral Kurtosis Heatmap for Interval')
    axs[1].set_ylabel('Time Blocks')
    axs[1].set_xlabel('Frequency (MHz)')

    # Add the colorbar showing actual SK values
    cbar = fig.colorbar(cax, ax=axs[1], orientation='horizontal', pad=0.1, aspect=40)

    # Add a legend for the color-coded segments
    axs[1].plot([], [], color='blue', lw=2, label='30% - 50% bad channels')
    axs[1].plot([], [], color='red', lw=2, label='>50% bad channels')
    axs[1].legend(loc='upper right', fontsize=8)

    # Save the combined plots
    output_file = os.path.join(output_folder, 'combined_sk_heatmap_and_histogram_interval.png')
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()

    # Save block data to a text file for this interval
    block_data_file = os.path.join(output_folder, 'block_bad_channel_percentages.txt')
    with open(block_data_file, 'w') as f:
        f.write('Block Number | Percentage of Bad Channels (%)\n')
        for block, percent in block_data:
            f.write(f'{block:8d} | {percent:6.2f}\n')

def run_analysis_for_intervals(file_path, output_folder, target_resolution_ms=1.0, num_intervals=1, interval_duration=100):
    sk_list = []
    sk_variance_list = []
    all_outliers_sk = []
    all_frequent_outliers_ranges = []

    # Get header information based on file type
    header = get_file_header(file_path)

    # Calculate observation times and coordinates
    start_mjd = header['STT_IMJD'] + (header['STT_SMJD'] + header['STT_OFFS']) / 86400.0
    observation_start_time = Time(start_mjd, format='mjd')
    
    total_duration = header['TBIN'] * header['NAXIS2']
    observation_end_time = observation_start_time + TimeDelta(total_duration, format='sec')

    # Calculate positions
    azimuth_start, elevation_start = calculate_az_el(header, observation_start_time)
    azimuth_end, elevation_end = calculate_az_el(header, observation_end_time)
    # process_interval(0, file_path, output_folder, 
    #             target_resolution_ms, interval_duration, num_intervals,
    #             azimuth_start, elevation_start
    #         )
    # Process intervals in parallel
    with ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(
                process_interval, i, file_path, output_folder, 
                target_resolution_ms, interval_duration, num_intervals,
                azimuth_start, elevation_start
            )
            for i in range(num_intervals)
        ]

        for future in futures:
            sk, sk_variance, outliers_sk, frequent_outliers_sk, frequencies = future.result()
            sk_list.append(sk)
            sk_variance_list.append(sk_variance)
            all_outliers_sk.append(outliers_sk)
            all_frequent_outliers_ranges.extend(frequent_outliers_sk)

    # Combine results
    concatenated_sk = np.concatenate(sk_list, axis=0)
    concatenated_sk_variance = np.concatenate(sk_variance_list, axis=0)
    concatenated_outliers_sk = np.concatenate(all_outliers_sk, axis=0)
    merged_frequent_outliers = merge_frequency_ranges(all_frequent_outliers_ranges)

    # Save final results
    combined_outliers_file = os.path.join(output_folder, 'combined_frequent_outliers.txt')
    save_outliers_to_file(merged_frequent_outliers, combined_outliers_file, is_freq=True)

    plot_combined_sk_heatmap_and_histogram(
        concatenated_sk, concatenated_sk_variance, frequencies,
        concatenated_outliers_sk, merged_frequent_outliers, output_folder,
        azimuth_start, elevation_start, azimuth_end, elevation_end,
        header, file_path
    )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process FITS or filterbank file and analyze RFI.')
    parser.add_argument('input_file', type=str, help='Path to the FITS or filterbank file.')
    parser.add_argument('output_folder', type=str, help='Output folder for plots and results.')
    parser.add_argument('--target_resolution_ms', type=float, default=1.0,
                      help='Target time resolution in milliseconds for downsampling.')
    parser.add_argument('--num_intervals', type=int, default=18,
                      help='Number of 100-second intervals to process.')

    args = parser.parse_args()
    run_analysis_for_intervals(args.input_file, args.output_folder,
                             args.target_resolution_ms, args.num_intervals)

