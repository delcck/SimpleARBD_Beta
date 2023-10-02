#! /usr/bin/env python2
# /software/python3-3.4.1/bin/python3
#! /usr/bin/env python3

from __future__ import absolute_import, print_function
import numpy as np
from scipy import signal
import sys
from os.path import isfile, getmtime
from gridData import Grid
import math
import os

def info(*obj):
    print('INFO: ',obj , file=sys.stderr)

def loadGrid(file, **kwargs):
    return Grid(file)

    ## Load cached grid where possible (ONLY use this if pickles are trusted)
    datafile = '{}.pickle'.format(file)
    assert(isfile(file))

    ## Read cache or original file
    if isfile(datafile) and getmtime(datafile) > getmtime(file):
        info("Reading {}".format(datafile))
        data = Grid()
        data.load(datafile)
    else:
        info("Creating {}".format(datafile))
        data = Grid(file)
        data.save(datafile)
    return data

def Get_damping_coefficients(hydroproFile, massFile, inertiaFile, outFile):

    lineNum = 1
    with open(hydroproFile, 'r') as fin:
        ## skip 49 lines
        while lineNum <= 48:
            fin.readline()
            lineNum += 1

        ## read 3 lines
        Dx = float(fin.readline().strip().split()[0])
        Dy = float(fin.readline().strip().split()[1])
        Dz = float(fin.readline().strip().split()[2])

        ## skip 2 lines
        fin.readline()
        fin.readline()

        ## read 3 lines
        Rx = float(fin.readline().strip().split()[3])
        Ry = float(fin.readline().strip().split()[4])
        Rz = float(fin.readline().strip().split()[5])

    with open(massFile ,'r') as fin:
        mass = float(fin.readline().strip())
    with open(inertiaFile, 'r') as fin:
        inertia = [float(x) for x in fin.readline().strip().split()]

    ## convert
    # units "(295 k K) / (( cm^2/s) *  amu)" "1/ns"
    Dx, Dy, Dz = [24.527692/(x*mass) for x in [Dx,Dy,Dz]]
    #print(Dx,Dy,Dz)

    # units "(295 k K) / ((1 /s) *  amu AA^2)" "1/ns"
    Rx, Ry, Rz = [2.4527692e+17 / (x*mass) for x, mass in zip([Rx,Ry,Rz],inertia)]
    #print(Rx,Ry,Rz)

    with open(outFile, 'w') as fout:
        fout.write(' '.join([str(Dx), str(Dy), str(Dz)]) + '\n')
        fout.write(' '.join([str(Rx), str(Ry), str(Rz)]))


def Fix_charge(inFile, outFile, netChargeFile):
    cmd_in = "sed -r 's/^([0-9]+)e/\1.0e/g; s/ ([0-9]+)e/ \1.0e/' " + inFile + " > fix_charge_temp0.dx"
    os.system(cmd_in)
    cmd_in = "sed -r 's/^(-[0-9]+)e/\1.0e/g; s/ (-[0-9]+)e/ \1.0e/' fix_charge_temp0.dx > fix_charge_temp1.dx"
    os.system(cmd_in)

    resolution = 2

    with open(netChargeFile, 'r') as fout:
        netCharge = float(fout.readline().strip())

    ## Load Data
    #g = loadGrid(inFile)
    g = loadGrid('fix_charge_temp1.dx')
    g.grid = g.grid * resolution**3

    ## Apply upper and lower bounds
    ids = np.where(np.abs(g.grid[:]) > 0.01)

    numPoints = np.size(ids)
    info(np.sum(g.grid), numPoints, np.sum(g.grid)/numPoints)

    ## Remove excess charge (in loop due to machine error)
    while np.abs(np.sum(g.grid) - netCharge) > 0.0001:
        g.grid[ids] = g.grid[ids] + (netCharge-np.sum(g.grid))/numPoints
        info(np.sum(g.grid), numPoints, np.sum(g.grid)/numPoints)

    info("Final charge", np.sum(g.grid))

    ## Write output
    g.export(outFile)

def Bound_grid(inFile, outFile, lowerBound, upperBound):
    cmd_in = "sed -r 's/^([0-9]+)e/\1.0e/g; s/ ([0-9]+)e/ \1.0e/' " + inFile + " > bound_grid_temp0.dx"
    os.system(cmd_in)
    cmd_in = "sed -r 's/^(-[0-9]+)e/\1.0e/g; s/ (-[0-9]+)e/ \1.0e/' bound_grid_temp0.dx > bound_grid_temp1.dx"
    os.system(cmd_in)

    assert(lowerBound < upperBound)

    ## Load Data
    #g = loadGrid(inFile)
    g = loadGrid('bound_grid_temp1.dx')

    ## Apply upper and lower bounds
    ids = np.where(g.grid > upperBound)
    g.grid[ids] = upperBound
    ids = np.where(g.grid < lowerBound)
    g.grid[ids] = lowerBound

    ## Write output
    g.export(outFile)

def blur3Dgrid(g,blur):
    sideLen = 2*int(blur*3)+1
    gauss = signal.gaussian( sideLen, blur )
    i = np.arange( sideLen )
    i,j,k = np.meshgrid(i,i,i)
    kernel = gauss[i]*gauss[j]*gauss[k]
    kernel = kernel/kernel.sum()
    return signal.fftconvolve(g,kernel,mode='same')

def Find_segments_num(dimensions, threshold=300):
    in_xyz = [float(elm) for elm in dimensions]
    segments = [math.ceil(elm / threshold) for elm in in_xyz]
    return segments[0], segments[1], segments[2]

def Find_boundary_resolution(cellBasisVector1=[10,0,0],
                             cellBasisVector2=[0,10,0],
                             cellBasisVector3=[0,0,10],
                             resolution=2):
    min_PBC_length = min([max(cellBasisVector1), max(cellBasisVector2), max(cellBasisVector3)])
    if min_PBC_length < resolution:
        dx = round(min_PBC_length / 2)
    else:
        dx = resolution

    n1 = round(np.linalg.norm(cellBasisVector1)/dx)
    n2 = round(np.linalg.norm(cellBasisVector2)/dx)
    n3 = round(np.linalg.norm(cellBasisVector3)/dx)

    return dx, n1, n2, n3

def Find_boundary_end_points(cellBasisVector1=[10,0,0],
                             cellBasisVector2=[0,10,0],
                             cellBasisVector3=[0,0,10],
                             cellOrigin=[0,0,0]):
    xyz_travel = [sum([cellBasisVector1[ind], cellBasisVector2[ind], cellBasisVector3[ind]]) for ind in range(3)]
    xyz_min = [cellOrigin[ind] - xyz_travel[ind] * 0.5 for ind in range(3)]
    xyz_max = [cellOrigin[ind] + xyz_travel[ind] * 0.5 for ind in range(3)]

    return (xyz_min[0], xyz_min[1], xyz_min[2], xyz_max[0], xyz_max[1], xyz_max[2])

def Create_a_rectangular_mesh(wallRangeX, wallRangeY, wallRangeZ, blur, dd):

    x = np.linspace(wallRangeX[0], wallRangeX[1], int((wallRangeX[1]-wallRangeX[0])/dd) + 1)
    y = np.linspace(wallRangeY[0], wallRangeY[1], int((wallRangeY[1]-wallRangeY[0])/dd) + 1)
    z = np.linspace(wallRangeZ[0], wallRangeZ[1], int((wallRangeZ[1]-wallRangeZ[0])/dd) + 1)
    dx = np.mean(np.diff(x))
    dy = np.mean(np.diff(y))
    dz = np.mean(np.diff(z))

    X,Y,Z = np.meshgrid(x,y,z,indexing='ij')

    return X, Y, Z, dx, dy, dz

def Create_boundary_region(mx, my, mz,
                           n1, n2, n3,
                           cellBasisVector1=[10,0,0],
                           cellBasisVector2=[0,10,0],
                           cellBasisVector3=[0,0,10]):
    v1 = np.array(cellBasisVector1) / n1
    v2 = np.array(cellBasisVector2) / n2
    v3 = np.array(cellBasisVector3) / n3
    region_pts = []
    start_pt = np.array([mx, my, mz])
    for i in range(n1 + 1):
        for j in range(n2 + 1):
            for k in range(n3 + 1):
                region_pts.append(start_pt + i * v1 + j * v2 + k * v3)

    return region_pts

def Convert_math_pt_to_mesh_pt(region_pts, min_wallX, min_wallY, min_wallZ, dx, dy, dz):
    mesh_pts = []
    for pt in region_pts:
        indx = round((pt[0] - min_wallX) / dx)
        indy = round((pt[1] - min_wallY) / dy)
        indz = round((pt[2] - min_wallZ) / dz)
        mesh_pts.append((indx, indy, indz))
    return mesh_pts

def Create_the_well(X, wellDepth, mesh_pts, blur, dx, dy, dz, wX, wY, wZ, out_path):

    pot = np.zeros(np.shape(X))
    for pt in mesh_pts:
        pot[pt[0], pt[1], pt[2]] = wellDepth

    dd = np.mean([dx, dy, dz])

    pot_blur = blur3Dgrid(pot, blur/dd)

    origin = [wX, wY, wZ]
    opts = dict(delta=dd, origin=origin)
    g = Grid(pot_blur, **opts)
    g.export(out_path)

def Generate_spanning_vectors(bv1, bv2, bv3, dimensions, buff=5):

    dd = max(dimensions) + 2 * buff

    n1 = round(np.linalg.norm(bv1)/dd) + 1
    n2 = round(np.linalg.norm(bv2)/dd) + 1
    n3 = round(np.linalg.norm(bv3)/dd) + 1


    v1 = np.array(bv1) /n1
    v2 = np.array(bv2) /n2
    v3 = np.array(bv3) /n3


    return v1, v2, v3, n1, n2, n3

def Generate_coordinates(bv1, bv2, bv3, n1, n2, n3, num_copy, origin, replica_index):

    ori_vec = np.array(origin)

    if n1 * n2 * n3 > num_copy:
        check = True
    else:
        check = False

    np.random.seed(42 + replica_index)

    count = 0
    inds = []
    while count < num_copy:
        ind = (np.random.randint(0, n1), np.random.randint(0, n2), np.random.randint(0, n3))
        if not check:
            inds.append(ind)
            count += 1
        elif check:
            if not ind in inds:
                inds.append(ind)
                count += 1
            else:
                while ind in inds:
                    ind = (np.random.randint(0, n1), np.random.randint(0, n2), np.random.randint(0, n3))
                inds.append(ind)
                count += 1

    coors = []
    for ind in inds:
        r1 = (ind[0] - (n1 - 1) * 0.5) * bv1
        r2 = (ind[1] - (n2 - 1) * 0.5) * bv2
        r3 = (ind[2] - (n3 - 1) * 0.5) * bv3
        coors.append((r1 + r2 + r3 + ori_vec).tolist())

    return coors

def Create_null(grid_path='null.dx'):
    ones = np.ones([2,2,2])
    zeros = np.zeros([2,2,2])
    opts = dict(delta=3000,origin=-1500*np.array((1,1,1)))
    g = Grid(zeros,**opts)
    g.export(grid_path)
