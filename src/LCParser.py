import numpy as np
from pathlib import Path 
from typing import Optional, Dict
from scipy.signal import savgol_filter, peak_prominences
from src.utils import *
import matplotlib.pyplot as plt

class Parser:
    def __init__(self, file_path):
        self.file_path = file_path
    """
    Load the files and returns the raw contents
    """
    def __parse__(self):
        try:
            with open(self.file_path, 'r') as f:
                data = f.readlines()
        except FileNotFoundError:
            print(f"Error: File not found at {self.file_path}")
        except IOError:
            print(f"Error: Unable to read file at {self.file_path}")
        return data


class LCParser(Parser):
    def __init__(self, file_path):
        super().__init__(file_path)

        self.metadata = {
            "injection": [],
            "chromatogram_data": [],
            "signal_parameter": [],
        }
        # This will be a np array with [time, step, value] for the columns
        self.raw_data = [] 
        self.times = None
        self.values = None
        self.smoothed_values = None

        self.parse()
   
    def parse(self):
        """
        Load and parse the file specificed by file_path

        """
        # calls base class
        data = self.__parse__()
        start = 0

        # looks for start and end of Injection, Chroma Data and Signal Parameters and
        # parses those metadata into a dict by providing the data[start:end] into a parser
        for i, line in enumerate(data):
            
            if line.startswith("Injection Information:"):
                start = i

            elif line.startswith("Chromatogram Data Information:"):
                self.metadata["injection"] = parse_injection_info(data[start+1:i])
                start = i
            
            elif line.startswith("Signal Parameter Information:"):
                self.metadata["chromatogram_data"] = parse_chromatogram_data_info(data[start+1:i])
                start = i 

            elif line.startswith("Chromatogram Data"):
                self.metadata["signal_parameter"] = parse_signal_parameter_info(data[start+1:i])
                start = i
                break

        raw_chromatogram_data = data[start+2:]    
        # parse raw chroma data into an array
        self.raw_data = parse_raw_chromatogram_data(raw_chromatogram_data) 
        self.times = self.raw_data[:,0]
        self.values = self.raw_data[:,2]
    
    def peaks(self, height=0, apply_smoothing=True):
        """
        Returns peaks of the chroma data
    
        Parameters:
            height (float): minimum threshold to be considered a peak, by default we discard negative peaks.
            apply_smoothing(bool): whether to apply smoothing to the chroma values to reduce noisy peaks.
                this smooths the signal with scipy's savgol_filter
        
        Returns:
            peak_data: tuple of (left_threshold, right_threshold, height, peak_indice) for each peak so peak_data 
                is of length N, where:
                - left_threshold is the left time index of the peak
                - right_threshold is the right time index of the peak
                - height is the center height of the peak
                - peak_indice is the index of the peak in the raw data
        """
        peak_data = (None, None, None, None)
        if apply_smoothing:
            self.smoothed_values = self.filter_peaks()
            peak_data = self.__findpeaks__(self.smoothed_values, height)
        else:
            peak_data = self.__findpeaks__(self.values, height)
        return peak_data
        
    
    def __findpeaks__(self, values, height=0):  
        """
        helper function to find peaks in the raw chroma data
        """
        # finds indices in the values array where values[t] > values[t-1] &&  values[t] > values[t+1]
        is_peak = np.r_[True, values[1:] > values[:-1]] & np.r_[values[:-1] > values[1:], True] & (values > height)
        peak_indices = np.where(is_peak)[0]

        # don't consider the bounds of the data
        if 0 in peak_indices:
            peak_indices = np.delete(peak_indices, 0)
        if len(values)-1 in peak_indices:
            peak_indices = np.delete(len(values)-1, 0)


        heights = np.empty(peak_indices.shape[0], dtype=np.float64)
        left_thresh = np.empty(peak_indices.shape[0], dtype=np.intp)
        right_thresh = np.empty(peak_indices.shape[0], dtype=np.intp)

        # find left and right thresholds of the peaks as well as the height
        for i, peak_i in enumerate(peak_indices):
            left = peak_i
            min_left = values[left]
            # while (left >= 0) and (values[left] <= values[peak_i]):
            while (left >= 0) and (values[left] <= min_left):
                if values[left] < min_left:
                    min_left = values[left]
                    left_thresh[i] = left 
                left -= 1

            right = peak_i
            min_right = values[right]
            # while (right <= len(values)-1) and (values[right] <= values[peak_i]):
            while (right <= len(values)-1) and (values[right] <= min_right):
                if values[right] < min_right:
                    min_right = values[right]
                    right_thresh[i] = right
                right += 1

            # heights[i] = values[peak_i]
            heights[i] = values[peak_i] - max(min_left, min_right)
        
        return left_thresh, right_thresh, heights, peak_indices

    
    def filter_peaks(self, window=11, polyorder=1):
        """
        Smooth peaks with a savitzsky-golay filter

        Parameters:
            window (int): window size to for filter.
            polyorder (int): polynomial order to fit samples.
        
        Returns:
            filtered_signal (numpy.array): smoothed version of the peak values
        """
        filtered_signal = savgol_filter(self.values, window, polyorder)
        return filtered_signal

    def elution_volume(self, peak_data):
        """
        Calculate elution volume by integrating each of the peaks provided

        Parameters:
            peak_data (tuple): of left_thresholds, right_thresholds, heights for 
                all peaks detected       
        Returns:
            elution_volume (float): result of the integration
        """
        elution_volume = 0.0
        # use unsmoothed values to get more accurate value
        for left, right, _ in zip(*peak_data):
            peak_area = 0.0
            # use riemann sums to integrate
            for idx, peak_height in enumerate(self.values[left:right]):
                peak_area += max(peak_height, 0) * (self.times[left+idx] - self.times[left+idx-1])
            elution_volume += peak_area
        return elution_volume

    def get_times(self):
        return self.get_times
    
    def get_values(self):
        return self.values
    
    def get_smoothed_values(self):
        return self.smoothed_values
    
    @staticmethod
    def plot_signal(values, smoothed=None, peak_indices=None):
        """
        Plots the peak values.
        
        Parameters:
            values (numpy.ndarray): The 1D signal of values.
            smoothed (numpy.ndarray, optional): 1D signal of values after applying a smoothing function.
        peak_indices (list, optional): A list of indices where indicating peaks in the signal.
        """
        fig, ax = plt.subplots(figsize=(12, 6))
        
        ax.plot(values, label='Original Signal', color='blue')
        
        # Draw smoothed signal
        if smoothed is not None: 
            ax.plot(smoothed, label='Smoothed Signal', color='orange')

        # Draw open circles at the peaks
        if peak_indices is not None:
            if smoothed is not None:
                ax.plot(peak_indices, smoothed[peak_indices], 'ro', markerfacecolor='white')
            else:
                ax.plot(peak_indices, values[peak_indices], 'ro', markerfacecolor='white')

        ax.set_xlabel('Index')
        ax.set_ylabel('Value')
        plt.show()



def main(file_path):
    lc_parser = LCParser(file_path)
    left_thresh, right_thresh, peaks, peaks_ind = lc_parser.peaks(apply_smoothing=True)   
    peak_data = (left_thresh, right_thresh, peaks)
    elution_volume = lc_parser.elution_volume(peak_data)
    print(elution_volume, len(peaks_ind))
    import ipdb; ipdb.set_trace()    

    values = lc_parser.get_values()
    smoothed_values = lc_parser.get_smoothed_values()
    LCParser.plot_signal(values, smoothed_values, peaks_ind)


def main():
    parser = ArgumentParser(description='Run LCParser on a directory of data')
    parser.add_argument('--data_path,-d', type=str, required=True, help='Path to the data to run the function on. Either a single file or a flat folder of files')
    parser.add_argument('--detect', action='store_true')
    parser.add_argument('--smoothing', action='store_true', default=False, help='Enable smoothing')
    parser.add_argument('--save-plots', action='store_true', default=False, help='Save plots of LC data with detected peaks')


    args = parser.parse_args()
    smoothing = args.smoothing
    data_path = args.data_path


if __name__ == "__main__":
    data_dir = "/Users/jdeguzman/workspace/lc-parser/data"
    file_list = list(Path(data_dir).rglob("*.txt"))
    main(file_list[7])
    # for fpath in file_list:
    #     main(fpath)
