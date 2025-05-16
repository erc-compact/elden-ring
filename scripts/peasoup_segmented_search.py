import os


def peasoup_segmented_search(args):
    """
    Start multiple instance of peasoup according to the segmented search strategy.
    
    Args:
        Number of segements to divide the search space into.
        Number of samples in total
        FFT size for each segment
        
    Returns:
        None
    """
    
    # Number of segments
    num_segments = args.num_segments
    
    # Number of samples in total
    num_samples = args.nsamples
    
    # FFT size for each segment
    fft_size = args.fft_size
    
    # Calculate the start and end samples for each segment
    start_end_samples = start_and_end_sample(num_samples, num_segments)
    
    # Loop over the number of segments
    for i, (start_sample, end_sample) in enumerate(start_end_samples):
        # Calculate the number of samples in the segment
        samples_in_segment = end_sample - start_sample
        
        # Calculate the number of bins in the segment
        fft_size = fft_size / (i+1) 
        
        # Start peasoup for the segment
        
    


def start_and_end_sample():
    """
    Calculate the start and end sample for each segment.
    
    Args:
        Number of samples in total
        Number of segments
    Returns:
        List of tuples containing start and end samples for each segment
    """
    
    # Number of samples in total
    num_samples = args.nsamples
    
    # Number of segments
    num_segments = args.num_segments
    
    # Calculate the number of samples in each segment
    samples_per_segment = num_samples // num_segments
    
    # Initialize the list of tuples to store start and end samples
    start_end_samples = []
    
    # Loop over the number of segments
    for i in range(num_segments):
        # Calculate the start sample for the segment
        start_sample = i * samples_per_segment
        
        # Calculate the end sample for the segment
        end_sample = start_sample + samples_per_segment
        
        # Append the tuple to the list
        start_end_samples.append((start_sample, end_sample))
        
    return start_end_samples