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

def get_damping(t, x, N=100, threshold=0.1):
    # 1. Downsample via interpolation
    cs = sp.interpolate.CubicSpline(t, x)
    Z = cs(t)
    dt = t[2] - t[1]

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
    lambdahat = lambdas[0:M]

    # 8. Compute damping and frequency
    sks = np.log(lambdahat)
    alphas = np.real(sks)/dt
    omegas = np.imag(sks)/dt
    print(alphas)
    print(omegas)

    # 9. Compute amplitude
    Z2 = np.zeros((L, L), dtype=complex)
    for i in np.arange(0,L):  # L=len(lambdas)
        Z2[:,i] = np.power(lambdas[i], np.arange(0,L), dtype=complex).T
    h = np.linalg.lstsq(Z2, Z[:L], rcond=None)
    h = h[0]
    print(np.abs(h[:M]))

    return np.max(alphas)

params = {
    #0.50: np.arange(125.0, 150.0, 5.0),  # mach_number: [dynamic_pressures]
    #0.68: np.arange(115.0, 140.0, 5.0),
    # 0.90: np.arange(85.0, 110.0, 5.0),
    0.96: np.arange(65.0, 90.0, 5.0),
    # 1.07: np.arange(75.0, 95.0, 5.0),
    1.14: np.arange(205.0, 230.0, 5.0),
    }
modes = [1, 2, 3, 4]

main_path = os.getcwd()
figures_path = f'{main_path:s}/figures'
os.makedirs(figures_path, exist_ok=True)

for mach, qs in params.items():

    speed_index_paths = gl.glob(f'{main_path:s}/results/mach-{mach:.2f}/speed_index-*')

    # get speed indexes for this mach number
    speed_indexes = []
    for speed_index_path in speed_index_paths:
        speed_indexes.append(float(speed_index_path.split('/')[-1].split('-')[-1]))

    # loop through modal coordinates (DOFs)
    for mode_ind, mode in enumerate(modes):

        fig, ax = plt.subplots()

        # loop through speed indexes
        for speed_index_ind, speed_index in enumerate(speed_indexes):

            speed_index_path = speed_index_paths[speed_index_ind]

            filename = f'{speed_index_path:s}/unsteady/aehist_body1_mode{mode:d}.dat'
            ts, gdisp, gvel = load_mode(filename)

            damping = get_damping(ts, gdisp)
            ax.plot(ts, gdisp, label=f'{speed_index:.2f}, {damping:f}')
            ax.set_ylabel(f'Modal Displacement')
            ax.legend(loc='best', title='Speed index')
            ax.set_xlabel('Time (s)')

        ax.grid()
        fig.tight_layout()
        fig.savefig(f'{figures_path:s}/modal_displacement_mach-{mach:.2f}-mode-{mode:d}.png')
        plt.close()
        print(f'Saved modal displacement plot for mach = {mach:.2f}, mode = {mode:d}')
