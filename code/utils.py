from PIL import Image

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

from specutils import SpectralRegion
from astropy.modeling import models
from astropy import units as u
from specutils.spectra import Spectrum1D
from specutils.fitting import estimate_line_parameters
from specutils.fitting import fit_lines
from specutils.analysis import equivalent_width
from scipy.integrate import trapz

def dispose_zeros(x, y, del_zeros=True, clip=[None, None]):

    if del_zeros == True:
        bad_idx = list(np.where(y==0)[0])
        idxs = list(range(len(y)))
        good_idxs = [i for i in idxs if i not in bad_idx]
    else:
        good_idxs = list(range(len(y)))
    
    x = x[good_idxs][clip[0]:clip[1]]
    y = y[good_idxs][clip[0]:clip[1]]

    return x, y

    
def get_continuum(x, y, windows, smooth_factor=0.1, 
                  return_maxs=False, envelope_std=None):
    """
    Fit a spline to selected maxima points within bins on a spectrum.
    
    Parameters:
    x: array-like
        The x-axis data points (e.g., wavelength)
    y: array-like
        The y-axis data points (e.g., intensity)
    windows: int or list of tuples
        The bin sizes or a list of (start, end) indices with custom windows.
    smooth_factor: float or tuple
        Smoothing factor(s) for the spline fit.
    return_maxs: bool
        Whether to return the maxima points along with the spline.
    envelope_std: float or None
        Standard deviation factor for excluding points from the spline fit.
        
    Returns:
    tuple: Spline fitted data points and optionally maxima points
    """
    
    # Initialize smoothing factors
    if isinstance(smooth_factor, float):
        s1 = smooth_factor
        s2 = smooth_factor
    else:
        s1, s2 = smooth_factor

    # Initialize array to store maxima points
    y_binmax = np.full(y.shape, np.nan)

    # Find maxima points in windows
    if isinstance(windows, int):
        for i in range(0, len(y), windows):
            end = min(i+windows, len(y))
            max_idx = i + np.argmax(y[i:end])
            y_binmax[max_idx] = y[max_idx]
    else:
        idx_start = 0
        for idx_end, win in windows:
            if idx_end is None:
                idx_end = len(y)
            for i in range(idx_start, idx_end, win):
                end = min(i+win, len(y))
                max_idx = i + np.argmax(y[i:end])
                y_binmax[max_idx] = y[max_idx]
            idx_start = idx_end

    # Filter out invalid indices
    valid_indices = ~np.isnan(y_binmax)
    x_valid = x[valid_indices]
    y_valid = y_binmax[valid_indices]
    
    # Fit the first spline
    spline = UnivariateSpline(x_valid, y_valid, s=s1)

    # Optionally, filter maxima points based on standard deviation
    if envelope_std is not None:
        # Initialize array for residuals
        res_arr = np.full(y.shape, np.nan)
        
        # Compute residuals
        residuals = y_valid - spline(x_valid)
        std = np.std(residuals)
        
        k = 0
        for i, idx_bool in enumerate(valid_indices):
            if idx_bool:
                if residuals[k] > 0:
                    res_arr[i] = residuals[k] 
                elif residuals[k] < 0:
                    if np.abs(residuals[k]) < envelope_std*std:
                        res_arr[i] = residuals[k]
                    else:
                        y_binmax[i] = np.nan
                k += 1

        # Update valid indices
        valid_indices = ~np.isnan(res_arr)
        x_valid = x[valid_indices]
        y_valid = y_binmax[valid_indices]
        
        # Fit the second spline
        spline = UnivariateSpline(x_valid, y_valid, s=s2)

    return (spline(x), y_binmax) if return_maxs else spline(x)

def get_equivalent_width(x, y, x_rang, return_fit=False, model='gaussian', fit=True):

    spectrum = Spectrum1D(flux=y*u.Jy, spectral_axis=x*u.AA)

    # if y has negative values, make them 0
    y[y<0] = 0

    if fit==True:

        if model == 'gaussian':
            params = estimate_line_parameters(spectrum, models.Gaussian1D())
            amplitude = params.amplitude.value
            mean = params.mean.value
            stddev = params.stddev.value
            g_init = models.Gaussian1D(amplitude=amplitude*u.Jy, mean=mean*u.AA, stddev=stddev*u.AA)
            g_fit = fit_lines(spectrum, g_init)
        elif model == 'voigt':
            params = estimate_line_parameters(spectrum, models.Voigt1D())
            x_0 = params.x_0.value
            amplitude_L = params.amplitude_L.value
            fwhm_L = params.fwhm_L.value
            fwhm_G = params.fwhm_G.value
            g_init = models.Voigt1D(x_0=x_0*u.AA, amplitude_L=amplitude_L*u.Jy, fwhm_L=fwhm_L*u.AA, fwhm_G=fwhm_G*u.AA)
            g_fit = fit_lines(spectrum, g_init)
        elif model == 'lorentz':
            params = estimate_line_parameters(spectrum, models.Lorentz1D())
            x_0 = params.x_0.value
            amplitude = params.amplitude.value
            fwhm = params.fwhm.value
            g_init = models.Lorentz1D(x_0=x_0*u.AA, amplitude=amplitude*u.Jy, fwhm=fwhm*u.AA)
            g_fit = fit_lines(spectrum, g_init)
        elif model == 'moffat':
            params = estimate_line_parameters(spectrum, models.Moffat1D())
            amplitude = params.amplitude.value
            x_0 = params.x_0.value
            gamma = params.gamma.value
            alpha = params.alpha.value
            g_init = models.Moffat1D(amplitude=amplitude*u.Jy, x_0=x_0*u.AA, gamma=gamma*u.AA, alpha=alpha)

        else:
            raise ValueError('Model not recognized.')
        
        y_fit = g_fit(x*u.AA)

    else:
        y_fit = y


    
    spectrum_fit = Spectrum1D(flux=y_fit*u.Jy, spectral_axis=x*u.AA)

    ew = np.float32(equivalent_width(spectrum_fit, 
                                     regions=SpectralRegion(x_rang[0]*u.AA, x_rang[1]*u.AA)))

    ew = trapz(y, x)

    if return_fit:
        return ew, np.array(spectrum_fit.flux)
    else:
        return ew



def dict_swap(d, i, j):
    keys = list(d.keys())
    values = list(d.values())
    values[i], values[j] = values[j], values[i]
    d = dict(zip(keys, values))
    return d

def sort_dic(dic, lims=[None, None]):
    dic =  dict(sorted(dic.items(), key=lambda item: item[1]))

    dic = dict(list(dic.items())[lims[0]:lims[1]])

    return dic

def make_collage(paths, savepath, size=1080):

    # get files in savefold

    # Open each image and resize them to the same size
    
    images = [Image.open(path).resize((size, size)) for path in paths]

    # Create a new blank image to hold the collage
    collage = Image.new('RGB', (size*2, size*2))

    # Paste the four images into the collage
    collage.paste(images[0], (0, 0))
    collage.paste(images[1], (size, 0))
    collage.paste(images[2], (0, size))
    collage.paste(images[3], (size, size))

    # Save the collage image
    collage.save(savepath)

ordered_stellar_spectra = [
    'HD190429A_O4If',
    'HD46223_O4V',
    'HD91824_O7V',
    'HD69464_O7Ib(f)',
    'ALS_18929_B0V',
    'HD150898_B0Ib',
    'HD13267_B5Ia',
    'HD3369_B5V',
    'HD31295_A0V',
    'HD87737_A0I',
    'HD17378_A5Ia',
    'HD23194_A5V',
    'HD37227_F0II',
    'HD8829_F0V',
    'HD20902_F5I',
    'HD27524_F5V',
    'HD141004_G0V',
    'HD16901_G0I',
    'HD20630_G5V',
    'HD6474_G4Ia',
    'HD3651_K0V',
    'HD12014_K0Ib',
    'HD219978_K4p5Ib',
    'HD201091_K5V',
    'HD132933_M0I',
    'HD79211_M0V',
    'HD175588_M4II',
    'BD+195116_M4p5V'
]

# revers elist above
ordered_stellar_spectra = ordered_stellar_spectra[::-1]

def color_spectra(name, cmap='jet', range=(0,1)):

    spectra = ['O', 'B', 'A', 'F', 'G', 'K', 'M']
    vals = np.linspace(range[0], range[1], len(spectra))


    cmap = plt.get_cmap(cmap)
    colors_type = {s: (*cmap(v)[:3], 1) for s, v in zip(spectra, vals)}

    spectra_fulltype = name.split('_')[-1] 
    spectra_type = spectra_fulltype[0]
    spectra_subtype = spectra_fulltype[1]
    return  colors_type[spectra_type]

# spectral_lines = {
#     r'H$\alpha$': 6563,  # Balmer series
#     r'H$\beta$': 4861,   # Balmer series
#     r'H$\gamma$': 4340,  # Balmer series
#     r'H$\delta$': 4102,  # Balmer series
#     'HeI (1)': 4471,
#     'HeI (2)': 5876,
#     'HeI (3)': 6678,
#     'HeI (4)': 7065,
#     'HeII (1)': 4686,
#     'HeII (2)': 5411,
#     'LiI': 6708,
#     'CaI': 4226,
#     'CaII (1)': 3934,
#     'CaII (2)': 3969,
#     'CaII (3)': 8498,
#     'CaII (4)': 8542,
#     'CaII (5)': 8662,
#     'FeI (1)': 4046,
#     'FeI (2)': 4144,
#     'FeI (3)': 4325,
#     'FeI (4)': 4384,
#     'FeI (5)': 4405,
#     'FeII (1)': 4172,
#     'FeII (2)': 4233,
#     'FeII (3)': 4924,
#     'FeII (4)': 5018,
#     'FeII (5)': 5169,
#     'MgI': 5173,
#     'MgII': 4481,
#     'NaI (1)': 5890,
#     'NaI (2)': 5896,
#     'KI (1)': 7665,
#     'KI (2)': 7699,
#     'SiII': 4128,
#     'SiIII': 4552,
#     'SiIV': 4089,
#     'TiO (1)': 6158,
#     'TiO (2)': 6188,
#     'TiO (3)': 6270,
#     'TiO (4)': 6651,
#     'TiO (5)': 6678,
#     'TiO (6)': 7054,
#     'CII': 4267,
#     'OII': 4072,
#     'SII (1)': 6716,
#     'SII (2)': 6731,
#     'NII (1)': 6548,
#     'NII (2)': 6583
# }

# balmer_lines = {
#     r'H$\alpha$': 6563,  # n=3 to n=2
#     r'H$\beta$': 4861,   # n=4 to n=2
#     r'H$\gamma$': 4340,  # n=5 to n=2
#     r'H$\delta$': 4102,  # n=6 to n=2
#     r'H$\epsilon$': 3970, # n=7 to n=2
#     r'H$\zeta$': 3889,   # n=8 to n=2
#     r'H$\eta$': 3835,    # n=9 to n=2
#     r'H$\theta$': 3798,  # n=10 to n=2
# }

# # Ordering the dictionary by increasing wavelength
# spectral_lines = dict(sorted(spectral_lines.items(), key=lambda item: item[1]))
