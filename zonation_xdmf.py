# SUMMARY:      zonation_xdmf.py
# USAGE:        Generate XDMF file for eSTOMP zonation file
# ORG:          Pacific Northwest National Laboratory
# AUTHOR:       Xuehang Song
# E-MAIL:       xuehang.song@pnnl.gov
# ORIG-DATE:    July-2018
# DESCRIPTION:
# DESCRIP-END.
# COMMENTS:     only deal cartesian sturtured grids
#
# Last Change: 2018-07-02


import numpy as np
import re
import math
import glob


def length_conversion(x):
    """
    convert length units to m
    """

    return {
        'a': 1e-10,
        'ang': 1e-10,
        'angstrom': 1e-10,
        'ao': 1e-10,
        'cm': 0.01,
        'ffl': 109.728,
        'ft': 0.3048,
        'furlong': 201.168,
        'm': 1,
        'mi': 1609.344,
        'mile': 1609.344,
        "mm": 0.001,
        'rod': 5.0292,
        'yd': 0.9144
    }.get(x, 1)


def retrieve_grids(grids_in):
    """
    read grid information from eSTOMP input
    """

    # read raw input deck
    input_file = open(grids_in, "r")
    estomp_input = input_file.readlines()
    input_file.close()

    # remove comments and blank lines in input deck
    estomp_input = [re.split('[#!\n]', x)[0] for x in estomp_input]
    estomp_input = [x for x in estomp_input if x]

    # locate start of grid card
    grid_line = [i for i, s in enumerate(estomp_input) if "~Grid" in s][0]

    if "cartesian" in estomp_input[grid_line + 1].lower():
        print("Cartesian grids")
    else:
        sys.exit("Unfortunately, this scripts only can deal with cartesian grids")

    # read nx, ny, nz
    nx, ny, nz = map(int, re.split('[,]', estomp_input[grid_line + 2])[0:3])
    grid_value = []
    iline = 3
    d_flag = 0
    # loop lines of etomp inputs until have enough entry for grids
    # while estomp_input[grid_line + iline][0] != "~":
    while len(grid_value) < (1 + nx + 1 + ny + 1 + nz):
        line_data = estomp_input[grid_line + iline].split(",")
        ndata = int(math.floor(len(line_data) / 2))
        for idata in range(ndata):
            if ("@" in line_data[idata * 2]):
                d_flag = 1
                temp_n, temp_d = line_data[idata * 2].split("@")
                grid_value += [float(temp_d) *
                               length_conversion(line_data[idata * 2 + 1])] * int(temp_n)
            else:
                grid_value += [float(line_data[idata * 2]) *
                               length_conversion(line_data[idata * 2 + 1])]
        iline += 1

    # assign flatten grids values to x, y, z
    if d_flag == 1:
        xo = grid_value[0]
        dx = np.asarray(grid_value[1:1 + nx])
        yo = grid_value[1 + nx]
        dy = np.asarray(grid_value[1 + nx + 1:1 + nx + 1 + ny])
        zo = grid_value[1 + nx + 1 + ny]
        dz = np.asarray(
            grid_value[1 + nx + 1 + ny + 1:1 + nx + 1 + ny + 1 + nz])
        x = xo + np.cumsum(dx) - 0.5 * dx
        y = yo + np.cumsum(dy) - 0.5 * dy
        z = zo + np.cumsum(dz) - 0.5 * dz
        xe = xo + sum(dx)
        ye = yo + sum(dy)
        ze = zo + sum(dz)
    else:
        xo = grid_value[0]
        xe = grid_value[nx]
        dx = np.diff(grid_value[0:(nx+1)])
        yo = grid_value[nx+1]
        ye = grid_value[nx+1+ny]
        dy = np.diff(grid_value[nx+1:(nx+1+ny+1)])
        zo = grid_value[nx+1+ny+1]
        ze = grid_value[nx+1+ny+1+nz]
        dz = np.diff(grid_value[nx+1+ny+1:(nx+1+ny+1+nz+1)])
        x = xo + np.cumsum(dx) - 0.5 * dx
        y = yo + np.cumsum(dy) - 0.5 * dy
        z = zo + np.cumsum(dz) - 0.5 * dz
    print("Grid retrived from eSTOMP input")
    return xo, yo, zo, xe, ye, ze, dx, dy, dz, nx, ny, nz, x, y, z


simu_dir = '/pic/scratch/song884/dust/fy2018/by_10a/l2/'
input_file = simu_dir + "input"

ijk_files = {}
ijk_files["sat1"] = simu_dir + "ijk_sat1.txt"
ijk_files["sat2"] = simu_dir + "ijk_sat2.txt"
ijk_files["sat3"] = simu_dir + "ijk_sat3.txt"
ijk_files["n"] = simu_dir + "ijk_n.txt"
ijk_files["k"] = simu_dir + "ijk_k.txt"
ijk_files["rho"] = simu_dir + "ijk_rho.txt"
ijk_files["zonation"] = simu_dir + "nonuniform-zonation-by"

# get model coords
xo, yo, zo, xe, ye, ze, dx, dy, dz, nx, ny, nz, x, y, z = retrieve_grids(
    input_file)


# # remove comments and blank lines in input deck
# estomp_input = [re.split('[# ! \n]', x)[0] for x in estomp_input]
# estomp_input = [x for x in estomp_input if x]

# grid_line = [i for i, s in enumerate(estomp_input) if "~Grid" in s][0]
