from numpy import *
from numpy.linalg import *

step = 5

# This function returns the frequency of the largest powered signal


def getf(x):
    fourier = fft.rfft(x)
    mag = abs(fourier)**2
    maxa = argmax(mag[1:]) + 1
    r = fft.fftfreq(len(x), d=step)
    f = r[maxa]
    return f


def find_signal_peaks_and_troughs(y):
    """ y (signal) and identifies peaks and troughs with a defined minimum amplitude"""
    differrences = diff(y)
    peak_idx = []
    trough_idx = []

    current_state = None
    for idx, val in enumerate(differrences[1:], 1):
        previous_state = current_state

        if val < 0:
            current_state = "negative"

        elif val > 0:
            current_state = "positive"

        else:
            continue

        # A peak is when a positive preceeds a negative
        if current_state == "negative" and previous_state == "positive":
            peak_idx.append(idx)

        elif current_state == "positive" and previous_state == "negative":
            trough_idx.append(idx)

    return(peak_idx, trough_idx)


def get_amplitudes(peak_idx, trough_idx, y):
    """
    Takes arrays of peak and trough index values, y data to return the amplitude between each
    peak and trough
    """
    amplitudes = []

    for idx, i in enumerate(peak_idx):
        try:
            amp = abs(y[peak_idx[idx]] - y[trough_idx[idx]])
            amplitudes.append(amp)

        except IndexError:
            break

    return amplitudes


def normalise_signal(signal_data):
    max_signal = max(signal_data)
    min_signal = min(signal_data)

    def norm_function(x): return (x - min_signal)/(max_signal - min_signal)

    return norm_function(signal_data)


def count_intersects(signal_data):
    signal_median = median(signal_data)
    intersection_count = 0

    # Iterates over signal data, retrieving linear
    for i in range(0, len(signal_data)-1):
        start_y = signal_data[i]
        end_y = signal_data[i + 1]
        if (start_y < signal_median) and (end_y > signal_median):
            intersection_count = intersection_count + 1

        elif (start_y > signal_median) and (end_y < signal_median):
            intersection_count = intersection_count + 1

    return intersection_count


def average_value_difference(y, peaks_idx):
    values = [y[i] for i in peaks_idx]

    return(mean(diff(values)))


def compare_final_amp_to_all(final_amp, all_amps):

    diff_amps = []
    for amp in all_amps:
        diff_amps.append(final_amp - amp)

    return(abs(mean(diff_amps)))


def check_isnan(data):
    for idx, i in enumerate(data.T):
        if isnan(data[:, idx][0]):
            return True

    return False


def check_below_one(data):
    for idx, i in enumerate(data.T):
        if min(data[:, idx]) < 1:
            return True

    return False


def distance(data1, data2, parameters, model):
    # data1 is simulated, and has shape npoints x beta
    # data2 is real
    # model is the model number
    target_amp = 1e12 * 0.1
    target_period_freq = 2400
    # get data for competitor

    if check_isnan(data1) is True:
        return[None, None, None]

    if check_below_one(data1) is True:
        return[None, None, None]

    distances = []

    # Each strain has three distances
    for d in data1.T:
        d = d[500:]
        d_peak_idx, d_trough_idx = find_signal_peaks_and_troughs(d)
        d_amps = get_amplitudes(d_peak_idx, d_trough_idx, d)

        if len(d_amps) == 0:
            return [None, None, None]

        d_threshold_amps_count = 0
        period_freq = 1/getf(d)

        # Count number of amplitudes above the target
        for amp in d_amps:
            if amp > target_amp:
                d_threshold_amps_count = d_threshold_amps_count + 1

        d_final_amp = d_amps[-1]

        target_num_peaks = target_period_freq / period_freq

        dist_1 = abs(d_threshold_amps_count - target_num_peaks)
        dist_2 = d_final_amp - target_amp

        if dist_1 <= 0.9:
            dist_1 = 0

        if dist_2 > 0:
            dist_2 = 0

        else:
            dist_2 = abs(dist_2)

        if (period_freq < target_period_freq):
            dist_3 = 0

        else:
            dist_3 = abs(period_freq - target_period_freq)

        distances.extend([dist_1, dist_2, dist_3])

    return distances
