#!/usr/bin/env python

import os
import numpy as np
import glob as gl
import tecplot as tp
import matplotlib.pyplot as plt
import scipy as sp

fontsize = 12
plt.rc('font', size=fontsize)
plt.rc('figure', figsize=(12, 8))
plt.rc('axes', labelsize='large')
plt.rc('lines', linewidth=3)

def load_mode(filename):

    # Load modal displacement and velocity
    tempdata = tp.data.load_tecplot(
        filename,
        read_data_option=tp.constant.ReadDataOption.ReplaceInActiveFrame,
        initial_plot_first_zone_only=True,
        variables=['time', 'gdisp', 'gvel'],
    )

    # Rearrange data
    ts = tempdata.zone(0).values('time').as_numpy_array()
    gdisp = tempdata.zone(0).values('gdisp').as_numpy_array()
    gvel = tempdata.zone(0).values('gvel').as_numpy_array()

    return ts, gdisp, gvel

def get_damping_logdec(t, x):

    prominence = np.max(x)/2.0
    peaks = sp.signal.find_peaks(x, prominence=prominence)[0]
    n = len(peaks) - 2
    if n < 4:
        return 888

    z1 = x[peaks[1]]
    z2 = x[peaks[-2]]
    t0 = t[peaks[0]]
    t1 = t[peaks[1]]
    
    # damping ratio
    delta = (1.0/n)*np.log(z1/z2)
    zeta = 1/np.sqrt(1+(2*np.pi/delta)**2)

    # frequency
    wn = 2*np.pi/(t1-t0)

    return -1*zeta*wn

def get_damping(t, x, N=100, threshold=0.1, rho=100.0):

    # 1. Downsample via interpolation
    N = max(N, int(len(t)/10))
    cs = sp.interpolate.CubicSpline(t, x)
    ts = np.linspace(t[0], t[-1], N)
    Z = cs(ts)
    dt = ts[1] - ts[0]

    # 2. Set pencil parameter
    L = int(N/2 - 1)

    # 3. Fill Hankel matrix Y with samples Z
    Y = np.zeros((N-L, L+1))
    for i in np.arange(0, N-L):
        for j in np.arange(0, L+1):
            Y[i,j] = Z[i+j]
    
    # 4. Singular value decomposition
    U, S, Vt = np.linalg.svd(Y)

    # 5. Choose model order M
    smax = np.max(S)
    Snew = []
    for scurr in S:
        if scurr < threshold * smax:
            break
        Snew.append(scurr)
    M = len(Snew)

    # 6. Compute A matrix
    V = Vt.T
    Vhat = V[:,:M]
    Vhat1 = Vhat[:L, :]
    Vhat2 = Vhat[1:L+1, :]
    Vhat1t_pinv = np.linalg.pinv(Vhat1.T)
    Vhat2t = Vhat2.T
    A = Vhat1t_pinv @ Vhat2t

    # 7. Eigendecomposition
    lambdas = np.linalg.eigvals(A)
    lambdahat = lambdas[:M]

    # 8. Compute damping and frequency
    sks = np.log(lambdahat)
    alphas = np.real(sks)/dt
    omegas = np.imag(sks)/dt
    alpha_ks = (1/rho) * np.log(np.sum(np.exp(rho*alphas)))

    # 9. Compute amplitude
    Z2 = np.zeros((L, L), dtype=complex)
    for i in np.arange(0,L):  # L=len(lambdas)
        Z2[:,i] = np.power(lambdas[i], np.arange(0,L), dtype=complex).T
    h = np.linalg.lstsq(Z2, Z[:L], rcond=None)
    h = h[0]

    return alpha_ks

# MAIN SCRIPT -----------------------------------------------------------------

machs = [0.50, 0.68, 0.90, 0.96, 1.07, 1.14]  # -, mach number
modes = [1]

main_path = os.getcwd()
figures_path = f'{main_path:s}/figures'
os.makedirs(figures_path, exist_ok=True)

# Loop through Mach numbers
for mach in machs:

    # get speed indexes for this mach number
    speed_index_paths = gl.glob(f'{main_path:s}/results/mach-{mach:.2f}/speed_index-*')
    speed_indexes = []
    for speed_index_path in speed_index_paths:
        speed_indexes.append(float(speed_index_path.split('/')[-1].split('-')[-1]))

    # loop through modal coordinates (DOFs)
    for mode_ind, mode in enumerate(modes):

        fig, ax = plt.subplots()

        dampings = []

        # loop through speed indexes
        for speed_index_ind, speed_index in enumerate(speed_indexes):
            print(f'mach={mach}, speed index={speed_index}, mode={mode}')

            speed_index_path = speed_index_paths[speed_index_ind]
            filename = f'{speed_index_path:s}/unsteady/aehist_body1_mode{mode:d}.dat'

            # Load data
            ts, gdisp, gvel = load_mode(filename)

            # Calculate damping
            print(f'ts.shape:\t{ts.shape}\ngdisp.shape:\t{gdisp.shape}')
            tend = min(3000,len(ts))
            damping = get_damping(ts[200:tend], gdisp[200:tend], threshold=0.5, rho=10.0)
            dampings.append(damping)
            # damping_ld = get_damping_logdec(ts[200:tend], gdisp[200:tend])

            # Plotting
            ax.plot(ts[:3000], gdisp[:3000], label=f'{speed_index:.3f}, {damping:.3f}')
            ax.set_ylabel(f'Modal Displacement')
            ax.legend(loc='best', title='Speed index, damping')
            ax.set_xlabel('Time (s)')

        # Calculate FSI (via interpolation)
        cs = sp.interpolate.CubicSpline(speed_indexes, dampings)
        speed_index_search = np.linspace(min(speed_indexes), max(speed_indexes), 100)
        # damping_interps = np.interp(speed_index_search, speed_indexes, dampings)
        damping_interps = cs(speed_index_search)
        flutter_array_ind = np.absolute(damping_interps).argmin()
        fsi = speed_index_search[flutter_array_ind]
        ax.set_title(f'FSI = {fsi:.3f}')

        ax.grid()
        fig.tight_layout()
        fig.savefig(f'{figures_path:s}/modal_displacement_mach-{mach:.2f}-mode-{mode:d}.png')
        plt.close()
        print(f'Saved modal displacement plot for mach = {mach:.2f}, mode = {mode:d}')
