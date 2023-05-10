import numpy as np
from scipy.signal import savgol_filter as savgol
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt




def point_correlation_1(func1, func2, i, delta):
    return (delta**2 + (func1[i+1] - func1[i]) * (func2[i+1] - func2[i])) / \
        ( np.sqrt(delta**2 + (func1[i+1] - func1[i])**2) * np.sqrt(delta**2 + (func2[i+1] - func2[i])**2) )


def point_correlation_2(func1, func2, i, delta):
    theta_1 = np.arctan((func1[i+1] - func1[i]) / delta)
    theta_2 = np.arctan((func2[i+1] - func2[i]) / delta)
    diff = np.abs(theta_2 - theta_1)
    corr = (1 / diff) - (1 / (2 * np.pi))
    return corr
    

def correlations_self(func1, func2, x_data, corrfunc= point_correlation_1):
    return [corrfunc(func1, func2, i, delta) for i, delta in zip(range(len(x_data[:-1])), np.diff(x_data))]


def correlations_numpy(func1, func2):
    assert len(func1) == len(func2), "Arrays of unequal length!"
    corrs = [np.corrcoef([func1[i-2:i+3], func2[i-2:i+3]])[0,1] for i in range(2, len(func1)-2)]
    corrs_zinf = list(map(lambda x: (-2 / (x-1)) - 1, corrs))
    return corrs_zinf


def max_shft_correlation(func1, func2, x_data):
    mean_correlations = np.zeros(len(x_data))
    for i in range (len(x_data)):
        shift_func1 = deque(func1)
        shift_func1.rotate(i)
        mean_correlations[i] = np.mean(correlations_self(shift_func1, func2, x_data))
    return np.max(mean_correlations)



#function to work on later
def intrp_clip(arr, low_lim, up_lim, start_y, end_y):

    if (not isinstance(arr, np.ndarray)):
        arr = np.array(arr)
        
    correct_mask = (arr > low_lim) & (arr < up_lim)
    corrected_arr = arr
    for i in range(len(arr)):
        if (not correct_mask[i]):
            if (len(np.where(correct_mask[:i] == True)) == 0):
                prev_y = start_y
                prev_x = 0
            else :
                prev_y = arr[:i][correct_mask[:i]][-1]
                prev_x = np.argwhere(correct_mask[:i])[-1]
            if (len(np.where(correct_mask[i:] == True)) == 0):
                next_y = end_y
                next_x = len(arr) - 1
            else :
                next_y = arr[i:][correct_mask[i:]][0]
                next_x = np.argwhere(correct_mask[i:])[0]
            slope = (next_y - prev_y) / (next_x - prev_x)
            corrected_arr[i] = prev_y + slope * (i - prev_x)
        else :
            corrected_arr[i] = arr[i]
    return corrected_arr
