import xml.etree.ElementTree as ET
import argparse

def extract_frequencies(xml_file):
    """
    Extracts periods from the XML file and calculates frequencies (1/period).

    Args:
        xml_file (str): The path to the XML file.
        fft_size (int): The size of the FFT to compute the frequency bin spacing.

    Returns:
        List of tuples: List containing frequencies in Hz and the corresponding 1/fft_size values.
    """
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Find the period information
    periods = []
    for i, candidate in enumerate(root.findall('.//candidate')):
        period = float(candidate.find('period').text)  # Assuming period tag exists
        frequencies = 1.0 / period  # Frequency is the inverse of the period
        periods.append(frequencies)

    tsamp = float(root.find('.//tsamp').text)
    nsamples = int(root.find('.//nsamples').text)
    
    # Calculate bin spacing upto 4 decimal places
    bin_spacing = float(5.0 / (tsamp * nsamples))  # Nyquist frequency is 1/(2 * tsamp)

    # Prepare the list of frequencies with the corresponding bin spacing
    freq_bin_pairs = [(freq, bin_spacing) for freq in periods]

    return freq_bin_pairs

def write_birdies_to_file(frequencies, default_birdies_file="None"):
    """
    Write the frequencies and their corresponding 1/fft_size values to a file.

    Args:
        frequencies (list): List of tuples with frequency and bin spacing.
        output_file (str): The file to write the frequencies to.
    """
    with open("birdies.txt", 'w') as f:
        # add default birdies from default_birdies.txt
        try:
            with open(default_birdies_file, 'r') as df:
                for line in df:
                    f.write(line)
        except FileNotFoundError:
            pass  # If the default file doesn't exist, just skip it
        
        for freq, spacing in frequencies:
            f.write(f"{freq:.6f} {spacing:.6f}\n")

def main(args):
    frequencies = extract_frequencies(args.xml_file)
    write_birdies_to_file(frequencies, args.default_birdies_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract periods from an XML file, calculate frequencies, and save them in a specific format.")
    parser.add_argument('--xml_file', type=str, help='Path to the XML file.')
    parser.add_argument('--default_birdies_file', type=str, default="None", help='Path to the default birdies file.')
    args = parser.parse_args()
    main(args)
