from argparse import ArgumentParser
from pathlib import Path
from src.LCParser import LCParser

def run(file_path, smooth, save_plot=False):
    lc_parser = LCParser(file_path)

    # print metrics
    left_thresh, right_thresh, peaks, peaks_ind = lc_parser.peaks(apply_smoothing=smooth)   
    peak_data = (left_thresh, right_thresh, peaks)
    elution_volume = lc_parser.elution_volume(peak_data)
    
    import ipdb; ipdb.set_trace()
    # save plots
    if save_plot:
        values = lc_parser.get_values()
        smoothed_values = lc_parser.get_smoothed_values()
        LCParser.plot_signal(values, smoothed_values, peaks_ind)


def main():
    parser = ArgumentParser(description='Run LCParser on a directory of data')
    parser.add_argument('--data-path', type=str, required=True, help='Path to the data to run the function on. Either a single file or a flat folder of files')
    parser.add_argument('--smooth', action='store_true', default=False, help='optional: Enable smoothing')
    parser.add_argument('--save-plot', action='store_true', default=False, help='optional: Save plots of LC data with detected peaks')

    args = parser.parse_args()
    print(f"Received args: {args}")

    smooth = args.smooth
    data_path = Path(args.data_path)
    save_plot = args.save_plot

    if data_path.is_file():
        print(f"Running function on file")
        print(f"Processing file: {data_path}")
        run(data_path, smooth, save_plot)

    elif data_path.is_dir():
        print(f"Running function on directory")
        # Iterate through files in the directory
        for file_path in data_path.glob('*'):
            print(f"Processing file: {file_path}")
            # run(file_path, smooth, save_plot)
    else:
        print(f"Invalid data_path: {data_path}")
        return


if __name__ == "__main__":
    main()