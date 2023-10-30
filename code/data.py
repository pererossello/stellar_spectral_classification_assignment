import numpy as np


balmer_lines = {
    r'H$\alpha$': 6563,  # n=3 to n=2
    r'H$\beta$': 4861,   # n=4 to n=2
    r'H$\gamma$': 4340,  # n=5 to n=2
    r'H$\delta$': 4102,  # n=6 to n=2
    r'H$\epsilon$': 3970, # n=7 to n=2
    r'H$\zeta$': 3889,   # n=8 to n=2
    r'H$\eta$': 3835,    # n=9 to n=2
    r'H$\theta$': 3798,  # n=10 to n=2
    r'H$\iota$': 3770,   # n=11 to n=2
    r'H$\kappa$': 3750,  # n=12 to n=2
    r'H$\lambda$': 3734, # n=13 to n=2
    r'H$\mu$': 3722,     # n=14 to n=2
    r'H$\nu$': 3712,     # n=15 to n=2
    r'H$\xi$': 3704,     # n=16 to n=2
    'H\u03BF': 3698, # n=17 to n=2
    # r'H$\pi$': 3692      # n=18 to n=2    
}

caii_lines = {
    'CaII H': 3933.66,
    'CaII K': 3968.47,
    'CaII (3)': 8498,
    'CaII (4)': 8542,
    'CaII (5)': 8662,
}

fei_lines = {
    'FeI (1)': 4046,
    'FeI (2)': 4144,
    'FeI (3)': 4325,
    'FeI (4)': 4384,
    'FeI (5)': 4405,
}

feii_lines = {
    'FeII (1)': 4172,
    'FeII (2)': 4233,
    'FeII (3)': 4924,
    'FeII (4)': 5018,
    'FeII (5)': 5169,
}

cii_lines = {
    'CII': 4267,
}


spectral_lines_p2 = {
    r'H$\alpha$': 6563,  # Balmer series
    r'H$\beta$': 4861,   # Balmer series
    r'H$\gamma$': 4340,  # Balmer series
    r'H$\delta$': 4102,  # Balmer series
    r'H$\epsilon$': 3970, # n=7 to n=2
    r'H$\zeta$': 3889,   # n=8 to n=2
    r'H$\eta$': 3835,    # n=9 to n=2
    r'H$\theta$': 3798,  # n=10 to n=2
    'CaII H': 3933.66,
    'CaII K': 3968.47,
    'FeI (2)': 4144,
    'HeI (1)': 4121,
    'HeI (2)': 4387,
    'HeI (3)': 4471,
    'HeI (4)': 5876,
    'HeI (5)': 6678,
    'HeI (6)': 7065,
    'FeII (3)': 4924,
    'FeII (4)': 5018,
    'MgII': 4481,
    'NaI (1)': 5890,
    'NaI (2)': 5896,
    'SiII': 4128,
    'SiII': 4128,
    'SiIII': 4552,
    'SiIV': 4089,
    'CII': 4267,
    'OII': 4416,
    'HeII': 4686
}

spectral_lines_p2_balmer = {
    r'H$\alpha$': 6563,  # n=3 to n=2
    r'H$\beta$': 4861,   # n=4 to n=2
    r'H$\gamma$': 4340,  # n=5 to n=2
    r'H$\delta$': 4102,  # n=6 to n=2
    r'H$\epsilon$': 3970, # n=7 to n=2
    r'H$\zeta$': 3889,   # n=8 to n=2
    r'H$\eta$': 3835,    # n=9 to n=2
    r'H$\theta$': 3798,  # n=10 to n=2
    r'H$\iota$': 3770,   # n=11 to n=2
    r'H$\kappa$': 3750,  # n=12 to n=2
    r'H$\lambda$': 3734, # n=13 to n=2
    r'H$\mu$': 3722,     # n=14 to n=2
    r'H$\nu$': 3712,     # n=15 to n=2
    r'H$\xi$': 3704,     # n=16 to n=2
    'H\u03BF': 3698, # n=17 to n=2
    'CaII H': 3933.66,
    'CaII K': 3968.47,
    'FeI (2)': 4144,
    'HeI (1)': 4121,
    'HeI (2)': 4387,
    'HeI (3)': 4471,
    'HeI (4)': 5876,
    'HeI (5)': 6678,
    'HeI (6)': 7065,
    'FeII (3)': 4924,
    'FeII (4)': 5018,
    'MgII': 4481,
    'NaI (1)': 5890,
    'NaI (2)': 5896,
    'SiII': 4128,
    'SiIII': 4552,
    'SiIV': 4089,
    'CII': 4267,
    'OII': 4416,
}






spectral_lines_p1 = {
    r'H$\alpha$': 6563,  # Balmer series
    r'H$\beta$': 4861,   # Balmer series
    r'H$\gamma$': 4340,  # Balmer series
    r'H$\delta$': 4102,  # Balmer series
    #'CaI': 4226,
    'CaII (1)': 3934,
    'CaII (2)': 3969,
    'FeI (1)': 4046,
    'FeI (2)': 4144,
    'FeI (3)': 4215,
    'FeI (4)': 4325,
    'FeI (5)': 4384,
    'FeI (6)': 4405,
    'CaI': 4226,
    # 'MgII': 4481,
    'NaI (1)': 5890,
    'NaI (2)': 5896,
    # 'KI (1)': 7665,
    # 'KI (2)': 7699,
    # 'OII': 4072,
    # 'CN': 3874,
    # 'CN (2)': 3884,
    # 'CH (1)': 4300,
    # 'CH A-X': 3870,
    # 'CH C-X': 3906,
    # r'TiO$\alpha$ (1)': 4954,  # α system
    # r'TiO$\alpha$ (2)': 5035,  # α system
    # r'TiO$\alpha$ (3)': 5167,  # α system
    r'TiO$\alpha$ (1)': 4761,  # α system
    r'TiO$\alpha$ (2)': 4954,  # α system
    r'TiO$\alpha$ (3)': 5167,  # α system
    r'TiO$\alpha$ (4)': 5448,  # α system
    r'TiO$\alpha$ (5)': 6158,  # α system11

    # r'TiO$\gamma$ (1)': 6158,  # γ system
    # r'TiO$\gamma$ (2)': 6180,  # γ system
    # r'TiO$\gamma$ (3)': 6193,  # γ system
    # r'TiO$\gamma$ (4)': 6259,  # γ system
    # r'TiO$\gamma$ (5)': 6651,  # γ system
    # r'TiO$\gamma$ (6)': 6717,  # γ system
    # r'TiO$\gamma$ (7)': 6745,  # γ system
    # r'TiO$\gamma$ (8)': 7054,  # γ system
    # r'TiO$\gamma$ (9)': 7088,  # γ system
    # r'TiO$\gamma$ (10)': 7126,  # γ system
    r'MgI (1)': 5167.3,  # b2 line
    r'MgI (2)': 5172.7,  # b3 line
    r'MgI (3)': 5183.6,   # b4 line
    r'SrII': 4077,   # b5 line
    r'CH G': 4300, 
    r'CN (1)': 3870, 
    r'CN (2)': 3885, 
    r'CN (3)': 4215,
    # r'O[III] 4959': 4959,  # Blue-green line
    # r'O[III] 5007': 5007, 


}

# then another possible suspect could be a blend of metal lines or other molecular bands. 
