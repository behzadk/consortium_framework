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

    norm_function = lambda x: (x - min_signal)/(max_signal - min_signal)

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

def distance(data1, data2, parameters, model):
    # data1 is simulated, and has shape npoints x beta
    # data2 is real
    # model is the model number

    target_freq = 0.02
    target_num_peaks = 15    # Target is 9 because we are starting measurement from t = 100
    target_peak_val_diff = 0
    target_amp = 1e12 * 0.1
    target_period_freq = 2400
    # get data for competitor
    dataA_1 = data1[:, 0]
    dataA_2 = data1[:, 1]
    dataB_1 = data1[:, 2]
    dataB_2 = data1[:, 3]
    dataN_C = data1[:, 4]
    dataN_X = data1[:, 5]
    dataS = data1[:, 6]

    if isnan(dataN_C[0]) or isnan(dataN_X[0]):
        return [None, None, None]

    if (
        min(dataA_1) < 0 or min(dataA_2) < 0 or
        min(dataB_1) < 0 or min(dataB_2) < 0 or
        min(dataN_C) < 1 or min(dataN_X) < 1 or
        min(dataS) < 0
    ):
        return[None, None, None]

    dataN_C = dataN_C[500:]

    dataN_C_peak_idx, dataN_C_trough_idx = find_signal_peaks_and_troughs(dataN_C)
    dataN_C_amps = get_amplitudes(dataN_C_peak_idx, dataN_C_trough_idx, dataN_C)
    average_amp = mean(dataN_C_amps)

    # if period_freq >= 2400:
    #     return[None, None, None]

    if len(dataN_C_amps) == 0:
        return[None, None, None]

    dataN_C_threshold_amps_count = 0
    period_freq = 1/getf(dataN_C)

    # Count number of amplitudes above the target
    for amp in dataN_C_amps:
        if amp > target_amp:
            dataN_C_threshold_amps_count = dataN_C_threshold_amps_count + 1

    dataN_C_final_amp = dataN_C_amps[-1]

    target_num_peaks = target_period_freq / period_freq

    d1 = abs(dataN_C_threshold_amps_count - target_num_peaks)
    if d1 <= 0.9:
        d1 = 0

    else:
        d1 = abs(dataN_C_threshold_amps_count - target_num_peaks)

    if dataN_C_final_amp > target_amp:
        d2 = 0

    else:
        d2 = abs(dataN_C_final_amp - target_amp)

    # Is this a good idea? Motivation is that we will have one distance to 
    # encourage standard deviation away from the middle value, while the other
    # encourages low deviation in peak values. Will it balance out nicely?

    if period_freq < target_period_freq:
        d3 = 0

    else:
        d3 = abs(period_freq - target_period_freq)

    # d4 = log10(std(dataN_C_amps))

    return [d1, d2, d3]
