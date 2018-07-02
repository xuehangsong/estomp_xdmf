# this scripts only deal sturtured grid for now
import h5py as h5
import numpy as np
import re
import sys
import math
import glob
from xml.etree import ElementTree as ET
from xml.dom import minidom
# from datetime import datetime, timedelta


# def collect_grids():
#     input_file = open(simu_dir + "input", "r")
#     print("hello world")


def prettify(element):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ET.tostring(element, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

# convert length units to m


def length_conversion(x):
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


def retrieve_grids(simu_dir):
    # read raw input deck
    input_file = open(simu_dir + "input", "r")
    estomp_input = input_file.readlines()
    input_file.close()

    # remove comments and blank lines in input deck
    estomp_input = [re.split('[# ! \n]', x)[0] for x in estomp_input]
    estomp_input = [x for x in estomp_input if x]

    # locate start of grid card
    grid_line = [i for i, s in enumerate(estomp_input) if "~Grid" in s][0]

    if "Cartesian" in estomp_input[grid_line + 1]:
        print("Cartesian grids")
    else:
        sys.exit("Unfortunately, this scripts only can deal with cartesian grids")

    # read nx, ny, nz
    nx, ny, nz = map(int, re.split('[, ]', estomp_input[grid_line + 2])[0:3])
    grid_value = []
    iline = 3
    cell_surface_flag = 0
    # loop lines of etomp inputs until have enough entry for grids
    # while estomp_input[grid_line + iline][0] != "~":    
    while len(grid_value) < (1 + nx + 1 + ny + 1 + nz):
        line_data = estomp_input[grid_line + iline].split(",")
        ndata = int(math.floor(len(line_data) / 2))
        for idata in range(ndata):
            if ("@" in line_data[idata * 2]):
                temp_n, temp_d = line_data[idata * 2].split("@")
                grid_value += [float(temp_d) *
                               length_conversion(line_data[idata * 2 + 1])] * int(temp_n)
            else:
                cell_surface_flag = 1                
                grid_value += [float(line_data[idata * 2]) *
                               length_conversion(line_data[idata * 2 + 1])]
        iline += 1

    # assign flatten grids values to x, y, z        
    if cell_surface_flag == 0:
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


# retrieve and ouput variables times from hdf5 files
def retrieve_variable_time(simu_dir, time_unit):
    all_h5 = np.sort(glob.glob(simu_dir + "plot*h5block"))
    plot_h5 = h5.File(all_h5[0], "r")
    varis = list(plot_h5['Step#0']['Block'])
    varis_with_units = [x.split("::")[0]
             for x in list(plot_h5.attrs.keys()) if "units" in x]
    units = {}
    for ivari in varis_with_units:
        print(ivari)
        units[ivari] = plot_h5.attrs[ivari + "::units"].decode("UTF-8")
    plot_h5.close()
    times = []
    for i_h5 in all_h5:
        plot_h5 = h5.File(i_h5, "r")
        times.append(float(plot_h5.attrs['time step'].decode("UTF-8")))
        plot_h5.close()
    if time_unit in ("s", "sec", "second"):
        times = np.asarray(times)
    elif time_unit in ("m", "min", "minute"):
        times = np.asarray(times) / 60
    elif time_unit in ("h", "hr", "hour"):
        times = np.asarray(times) / 60 / 60
    elif time_unit in ("d", "day"):
        times = np.asarray(times) / 60 / 60 / 24
    elif time_unit in ("y", "yr", "year"):
        times = np.asarray(times) / 60 / 60 / 24 / 365.25
    else:
        print("ERROR: unknow time units!")
    ntime = len(times)
    print("Varibles and times retrived from eSTOMP plot files")
    return ntime, times, varis, units, all_h5


# if __name__ == '__main__':
simu_dir = "/Users/song884/Dropbox/DVZ/WMAC/test_paraview/s7_1a/2018/"
simu_dir = "/home/xhsong/Dropbox/DVZ/WMAC/test_paraview/s7_1a/2018/"
simu_dir = '/pic/scratch/song884/bcomplex/model/'
time_unit = "yr"
xo, yo, zo, xe, ye, ze, dx, dy, dz, nx, ny, nz, x, y, z = retrieve_grids(
    simu_dir)
ntime, times, varis, units, all_h5 = retrieve_variable_time(
    simu_dir, time_unit)
all_xdmf = [x.replace("h5block", "xdmf") for x in all_h5]

xml_root = ET.Element("Xdmf", Version="3.0")
xml_domain = ET.SubElement(xml_root, "Domain")

# geometry is defined by VXVYVZ mode, three vectors
# assume toplogy geometry is time invariant
xml_toplogoy = ET.SubElement(xml_domain, "Topology",
                             {'TopologyType': '3DRECTMesh',
                              'Dimensions': "{0} {1} {2}".format(nz, ny, nx)})
xml_geometry = ET.SubElement(xml_domain, 'Geometry',
                             {'GeometryType': "VXVYVZ"})
xml_geometry_x = ET.SubElement(xml_geometry, 'DataItem',
                               {'Dimensions': str(nx),
                                "NumberType": "Float",
                                "Precision": "4",
                                "Format": "XML"})
xml_geometry_x.text = np.array_str(x).strip("[]").replace("\n", " ")
xml_geometry_y = ET.SubElement(xml_geometry, 'DataItem',
                               {'Dimensions': str(ny),
                                "NumberType": "Float",
                                "Precision": "4",
                                "Format": "XML"})
xml_geometry_y.text = np.array_str(y).strip("[]").replace("\n", " ")
xml_geometry_z = ET.SubElement(xml_geometry, 'DataItem',
                               {'Dimensions': str(nz),
                                "NumberType": "Float",
                                "Precision": "4",
                                "Format": "XML"})
xml_geometry_z.text = np.array_str(z).strip("[]").replace("\n", " ")


time_grid = ET.SubElement(xml_domain, 'Grid',
                          {'Name': 'TimeSeries',
                           'GridType': 'Collection',
                           'CollectionType': 'Temporal'})

for itime in range(ntime):
    print(itime)
    grid = ET.SubElement(time_grid, "Grid",
                         {'Name': all_h5[itime].split("/")[-1],
                          'GridType': 'Uniform'})
    time = ET.SubElement(grid, "Time", {'Value': str(times[itime]),
                                        "TimeType": "Single"})
    # use same toplogy and geometry for all snapshots
    topology_ref = ET.SubElement(
        grid, "Topology", {"Reference": "/Xdmf/Domain/Topology"})
    geometry_ref = ET.SubElement(
        grid, "Geometry", {"Reference": "/Xdmf/Domain/Geometry"})
    scalars = [x for x in varis if "Velocity" not in x]
    vectors = [x for x in varis if "Velocity" in x]
    nscalar = len(scalars)
    nvector = len(vectors)
    if nscalar > 0:
        scalar_handle = {}
        scalar_data_handle = {}
        for iscalar in range(nscalar):
            scalar_handle[iscalar] = ET.SubElement(grid, "Attribute",
                                                   {"Name": scalars[iscalar],
                                                    "AttributeType": "Scalar",
                                                    "Center": "Node"})
            scalar_data_handle[iscalar] = ET.SubElement(scalar_handle[iscalar], "DataItem",
                                                        {"Format": "HDF",
                                                         "NumberType": "Float",
                                                         "Precision": "4",
                                                         "Dimensions": "{0} {1} {2}".format(nz, ny, nx)})
            scalar_data_handle[iscalar].text = all_h5[itime].split("/")[-1] + \
                ":/Step#0/Block/" + scalars[iscalar] + "/0"
    # currently assume only velocity is a vector
    if nvector > 0:
        velocity = ET.SubElement(grid, "Attribute",
                                 {"Name": "Velocity",
                                  "AttributeType": "Vector",
                                  "Center": "Node"})
        velocity_data = ET.SubElement(velocity, "DataItem",
                                      {"ItemType": "Function",
                                       "Function": "JOIN($0, $1, $2)",
                                       "NumberType": "Float",
                                       "Precision": "4",
                                       "Dimensions": "{0} {1} {2} {3}".format(nz, ny, nx, 3)})
        velocity_x = ET.SubElement(velocity, "DataItem",
                                   {"Format": "HDF",
                                    "NumberType": "Float",
                                    "Precision": "4",
                                    "Dimensions": "{0} {1} {2}".format(nz, ny, nx)})
        velocity_x.text = all_h5[itime].split("/")[-1] + \
            ":/Step#0/Block/" + \
            'X-Dir. Aqueous Darcy Velocity (Node Centered)' + "/0"
        velocity_y = ET.SubElement(velocity, "DataItem",
                                   {"Format": "HDF",
                                    "NumberType": "Float",
                                    "Precision": "4",
                                    "Dimensions": "{0} {1} {2}".format(nz, ny, nx)})
        velocity_y.text = all_h5[itime].split("/")[-1] + \
            ":/Step#0/Block/" + \
            'Y-Dir. Aqueous Darcy Velocity (Node Centered)' + "/0"
        velocity_z = ET.SubElement(velocity, "DataItem",
                                   {"Format": "HDF",
                                    "NumberType": "Float",
                                    "Precision": "4",
                                    "Dimensions": "{0} {1} {2}".format(nz, ny, nx)})
        velocity_z.text = all_h5[itime].split("/")[-1] + \
            ":/Step#0/Block/" + \
            'Z-Dir. Aqueous Darcy Velocity (Node Centered)' + "/0"

fname = simu_dir + "estomp.xdmf"
with open(fname, 'w') as f:
    f.write(prettify(xml_root))
