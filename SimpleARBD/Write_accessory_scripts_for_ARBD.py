#--This is a package for generating accessory scripts that the SimpleARBD package will source to compute simulations input
import math

def Write_align_tcl(out_path='1-align_Chun.tcl'):
    text = '''lassign $argv prefix
set seltext all

proc rotationIsRightHanded {R {tol 0.01}} {
    set x [coordtrans $R {1 0 0}]
    set y [coordtrans $R {0 1 0}]
    set z [coordtrans $R {0 0 1}]

    set l [veclength [vecsub $z [veccross $x $y]]]
    return [expr {$l < $tol}]
}

set ID [mol new $prefix.psf]
mol addfile $prefix.pdb waitfor all
set all [atomselect $ID all]
set sel [atomselect $ID $seltext]

## Center system on $sel
$all moveby [vecinvert [measure center $sel weight mass]]

set continue 1
while { $continue } {
    ## Get current moment of inertia to determine rotation to align
    lassign [measure inertia $sel moments] com principleAxes

    ## Convert 3x3 rotation to 4x4 vmd transformation
    set R [trans_from_rotate $principleAxes]
    ## Fix left-handed principle axes sometimes returned by 'measure inertia'
    if { ! [rotationIsRightHanded $R] } {
        puts "This was true"
        # puts "rotation $R is not right handed! Fixing!"
        set R [transmult {{1 0 0 0} {0 1 0 0} {0 0 -1 0} {0 0 0 1}} $R]
    }

    puts "My rotation is here: $R"
    puts "My second line is here: [lassign [measure inertia $sel moments] com principleAxes]"

    ## Apply rotation and check that it worked
    $all move $R

    ## Get current moment of inertia to determine rotation to align

    lassign [measure inertia $sel moments] com principleAxes moments
    puts $principleAxes
    set goodcount 0
    foreach x0 {{1 0 0} {0 1 0} {0 0 1}} {
        set x [coordtrans [trans_from_rotate $principleAxes] $x0]
        if {[veclength [vecsub $x $x0]] < 0.01} {
            incr goodcount
        }
    }
    if { $goodcount == 3 } {
        set continue 0
    }
}

## Write transformation matrix to return to original conformation
set ch [open $prefix.rotate-back.txt w]
foreach line [trans_to_rotate [transtranspose $R]] {
    puts $ch $line
}
close $ch

## Write out moments of inertia
set ms ""
foreach m $moments { lappend ms [veclength $m] }
set ch [open $prefix.inertia.txt w]
puts $ch $ms
close $ch

## Write out mass
set ch [open $prefix.mass.txt w]
puts $ch [measure sumweights $sel weight mass]
close $ch

## Write out psf, pdb of transformed selection
$sel writepdb $prefix.aligned.pdb
$sel writepsf $prefix.aligned.psf'''
    with open(out_path, 'w') as fout:
        fout.write(text)

def Write_hydropro_config(system, mass,
                          temperature=22.0,
                          viscosity=0.01,
                          solvent_density=1.0):
    #--temperature is in degree celsius
    #--viscosity in poises
    #--solvent density is in g/cm^3
    out_path = 'hydropro.dat'
    text = '''hydro
''' + system + '''.hydro
hydro.pdb
1
2.9,
6,
1.2,
3.0,
''' + str(temperature) + ''',
''' + str(viscosity) + ''',
''' + str(mass) + ''',
1.0,
''' + str(solvent_density) + '''
-1,
-1,
0,
1
*'''
    with open(out_path, 'w') as fout:
        fout.write(text)


def Write_charge_density_tcl(
    resolution=2,
    out_path="2-charge-density_Chun.tcl"
):
    text = '''lassign $argv prefix
set resolution ''' + str(resolution) + '''
set ID [mol new $prefix.psf]
mol addfile $prefix.pdb
set all [atomselect $ID all]
set netCharge [measure sumweights $all weight charge]

## Write out charge density
volmap density $all -o $prefix.chargeDensity.dx -res $resolution -weight charge
## Write out pqr for subsequent EM calculation
$all writepqr $prefix.pqr

set ch [open $prefix.netCharge.dat w]
puts $ch $netCharge
close $ch

set ch [open $prefix.dimension.dat w]
set minmax [measure minmax $all]
set x_dim [expr [lindex [lindex $minmax 1] 0] - [lindex [lindex $minmax 0] 0]]
set y_dim [expr [lindex [lindex $minmax 1] 1] - [lindex [lindex $minmax 0] 1]]
set z_dim [expr [lindex [lindex $minmax 1] 2] - [lindex [lindex $minmax 0] 2]]
puts $ch $x_dim
puts $ch $y_dim
puts $ch $z_dim
close $ch'''
    with open(out_path, 'w') as fout:
        fout.write(text)

def Write_apbs_config(
    system='test', in_xyz_dim=[100, 100, 100], conc=0.15, temperature=300, buffer=50,
    Large_system='Off'
):
    out_path = str(system) + '.apbs'
    xyz_dim=[str(elm) for elm in in_xyz_dim]
    buffer = buffer
    xyz_cg = [str(math.ceil(elm + buffer)) for elm in in_xyz_dim]
    if Large_system == 'Off':
        xyz_dime = xyz_cg
        center = 'mol 1'
    elif Large_system == 'On':
        #dividend = max([math.ceil(elm/300) for elm in in_xyz_dim])
        dividend = 2
        xyz_dime = [str(math.ceil((elm + buffer) / dividend)) for elm in in_xyz_dim]
        #center = '0 0 0'
        center = 'mol 1'
    xyz_fg = xyz_cg
    text = \
    '''read
mol pqr ''' + str(system) + '''.pqr
end
elec
mg-auto
dime ''' + ' '.join(xyz_dime) + '''
cglen ''' + ' '.join(xyz_cg) + '''
cgcent ''' + center + '''
fglen ''' + ' '.join(xyz_fg) + '''
fgcent ''' + center + '''
mol 1
npbe
bcfl sdh
srfm smol
chgm spl2
ion 1 ''' + str(conc) + ''' 2.0
ion -1 ''' + str(conc) + ''' 2.0
pdie  12.0
sdie  78.54
sdens  10.0
srad  1.4
swin  0.3
temp  ''' + str(temperature) + '''
gamma  0.105
calcenergy no
calcforce no
write pot dx ''' +str(system) + '''.elec.tmp
end
quit'''
    with open(out_path, 'w') as fout:
        fout.write(text)

def Write_vdw_den_pot_tcl(potResolution=1,
                          denResolution=2,
                          num_heavy_cluster=3,
                          python_path='/usr/bin/python3',
                          cluster_result_path='vdw_assignment.dat',
                          out_path="3-vdw_Chun.tcl",
                          parameter_folder='parameters',
                          diffuse_or_static='diffuse'
    ):
    if diffuse_or_static == 'diffuse':
        text = '''set prefixes $argv
set potResolution ''' + str(potResolution) + '''
set denResolution ''' + str(denResolution) + '''
set IDs ""
set minRadius 0.5

package require ilstools

ILStools::readcharmmparams [glob ''' + parameter_folder + '''/*]

set ljParms ""
set lj_hyd ""

foreach prefix $prefixes {
    set ID [mol new $prefix.psf]
    mol addfile $prefix.pdb
    lappend IDs $ID

    ILStools::assigncharmmparams $ID; # sets radius and occupancy to rmin and eps

    set all [atomselect $ID "noh"]
    set hyd [atomselect $ID "hydrogen"]
    append ljParms " [lsort -unique [$all get {radius occupancy}]]"
    append lj_hyd " [lsort -unique [$hyd get {radius occupancy}]]"
}

set ljParms [lsort -unique $ljParms]
set lj_hyd [lsort -unique $lj_hyd]
set ljTypes ""

set ch [open tmp.dat w]
foreach vals $ljParms {
    lassign $vals r e
    if {$r < $minRadius} {continue};

    set count 0
    set types ""
    foreach ID $IDs {
        set sel [atomselect $ID "noh and radius $r and occupancy \\"$e\\""]
        incr count [$sel num]
        append types " [lsort -unique [$sel get type]]"
    }
    set types [lsort -unique $types]
    lappend ljTypes $types
    puts $ch "$r $e $count"
}
close $ch

set ch [open hyd.dat w]
foreach vals $lj_hyd {
    lassign $vals r e
    if {$r < $minRadius} {continue};

    set count 0
    set types ""
    foreach ID $IDs {
        set sel [atomselect $ID "hydrogen and radius $r and occupancy \\"$e\\""]
        incr count [$sel num]
        append types " [lsort -unique [$sel get type]]"
    }
    set types [lsort -unique $types]
    lappend ljTypes $types
    puts $ch "$r $e $count"
}
close $ch
puts "My LJ types are: $ljTypes"

#########################################
## Run python to cluster LJ parameters ##
#########################################
set tmpFile tmp.py
set ch [open $tmpFile w]
puts $ch {
import numpy as np
from scipy.cluster import vq

numClusters = ''' + str(num_heavy_cluster) + '''

d = np.loadtxt('tmp.dat')
d_hyd = np.loadtxt('hyd.dat')

## build new dataset with 'count' (d[:,2]) entries of each value
d2 = [np.outer( np.ones((1,int(d[i,2]))) , d[i,:2] ) for i in range(d.shape[0])]
d2_hyd = [np.outer( np.ones((1,int(d_hyd[i,2]))) , d_hyd[i,:2] ) for i in range(d_hyd.shape[0])]

d2 = np.vstack( d2 )
d2w = vq.whiten( d2 )              # normalize features

d2_hyd = np.vstack( d2_hyd )
d2w_hyd = vq.whiten( d2_hyd )              # normalize features


ind = 0;
scalebase = d2w[ind,:];
while np.any(scalebase == 0):
  ind = ind + 1
  scalebase = d2w[ind,:]

scale = d2[ind,:] / d2w[ind,:]

ind = 0;
scalebase = d2w_hyd[ind,:];
while np.any(scalebase == 0):
  ind = ind + 1
  scalebase = d2w_hyd[ind,:]

scale_hyd = d2_hyd[ind,:] / d2w_hyd[ind,:]

## perform cluster analysis
np.random.seed(seed=42)
codeBook,dist = vq.kmeans(d2w , numClusters)
assignments, dists = vq.vq(d[:,:2], scale * codeBook)

codeBook_hyd,dist_hyd = vq.kmeans(d2w_hyd , 1)
assignments_hyd, dists_hyd = vq.vq(d_hyd[:,:2], scale_hyd * codeBook_hyd)

assignments_total = np.concatenate((assignments, assignments_hyd + numClusters) , axis=None)
codeBook_total = np.concatenate((scale*codeBook, scale_hyd*codeBook_hyd), axis=0)

print(" ".join(["%d" % a for a in assignments_total]))
print(" ".join(["%.3f" % c[0] for c in codeBook_total]))
print(" ".join(["%.3f" % c[1] for c in codeBook_total]))
}
close $ch

unset env(PYTHONHOME)
set ch [open "| ''' + python_path + ''' $tmpFile" r]
#set ch [open "|/Scr/cmaffeo2/anaconda3/bin/python $tmpFile" r]

gets $ch assignments; list
gets $ch newR; list
gets $ch newE; list

close $ch


## Find atom types that map to each LJ parameter cluster
for {set i 0} {$i < [llength $newR]} {incr i} {
    set typeArray($i) ""
}
foreach i $assignments t $ljTypes {
    append typeArray($i) " $t"
}

set ch [open ''' + cluster_result_path + ''' w]
set i 0
foreach r $newR e $newE {
    puts $ch "$r $e $typeArray($i)"
    incr i
}
close $ch

########################################
## Write out densities and potentials ##
########################################
## Loop over molecules
foreach prefix $prefixes ID $IDs {
    set all [atomselect $ID all]
    set minmax [measure minmax $all]

    lassign $minmax min max
    set min [vecsub $min {12 12 12}]
    set max [vecadd $max {12 12 12}]
    set minmaxPot "{$min} {$max}"

    ## loop over new LJ params
    set i 0
    foreach r $newR e $newE {
        puts "My r is: $r"
        puts "My e is: $e"
        ## write out density grid
        set sel [atomselect $ID "type $typeArray($i)"]
        volmap interp $sel -o $prefix.vdw$i.den.dx -res $denResolution
        ## write potential grid
        volmap ils $ID $minmaxPot -cutoff 12.0 -o $prefix.vdw$i.pot.dx -res $potResolution -subres 3 -probecoor {{0.01 0.01 0.01}} -probevdw "{$e $r}" -maxenergy 20 -orient 1
    incr i
    }
}
exit'''
    elif diffuse_or_static == 'static':
        text = '''set prefixes $argv

set potResolution ''' + str(potResolution) + '''
set ch [open ''' + cluster_result_path + ''' r]
set i 0
set newR ""
set newE ""
while {[expr ![eof $ch]]} {
    gets $ch inputData($i)
    if {![string equal $inputData($i) ""]} {
        puts $inputData($i)
        lappend newR [lindex $inputData($i) 0]
        lappend newE [lindex $inputData($i) 1]
        incr i
    }
}
close $ch

package require ilstools
ILStools::readcharmmparams [glob ''' + parameter_folder + '''/*]

foreach prefix $prefixes {
    set ID [mol new $prefix.psf]
    mol addfile $prefix.pdb
    set all [atomselect top all]
    set minmax [measure minmax $all]

    lassign $minmax min max
    set min [vecsub $min {12 12 12}]
    set max [vecadd $max {12 12 12}]
    set minmaxPot "{$min} {$max}"

    ILStools::assigncharmmparams $ID;

    set i 0
    foreach r $newR e $newE {
        puts "My r is: $r"
        puts "My e is: $e"
        ## write potential grid
        volmap ils $ID $minmaxPot -cutoff 12.0 -o $prefix.vdw$i.pot.dx -res $potResolution -subres 3 -probecoor {{0.01 0.01 0.01}} -probevdw "{$e $r}" -maxenergy 20 -orient 1 -first 0 -last 0
        incr i
    }
}
exit'''
    with open(out_path, 'w') as fout:
        fout.write(text)

def Write_smoothing_tcl(in_file, out_file, GaussianWidth=2.5,
                       out_path='4-smooth_Chun.tcl'):
    text = \
    '''voltool smooth -sigma ''' + str(GaussianWidth) +\
    ''' -i ''' + str(in_file) + ''' -o ''' + str(out_file)
    with open(out_path, 'w') as fout:
        fout.write(text)

def Write_BD_configs(timestep='200e-6',
                     steps=20000000,
                     temperature=300,
                     num_heavy_cluster=3,
                     interactive='#',
                     grid_path='null.dx',
                     diffusible_objects=['test'],
                     num_copies_per_objects=[25],
                     masses=[100],
                     inertias=[[101,102,103]],
                     transDampings=[[104,105,106]],
                     rotDampings=[[107,108,109]],
                     static_objects=['static'],
                     Extra_pots_tags=['boundary.dx, vdw2'],
                     BDCoordinates=['bdcoord.txt'],
                     out_path='test.bd'):
    text = '''# seed 993
timestep ''' + str(timestep) + '''
steps ''' + str(steps) + '''
numberFluct 0
interparticleForce 1
tabulatedPotential 1
fullLongRange 0
temperature ''' + str(temperature) + '''
outputPeriod 1000
outputEnergyPeriod 1000
# make decompPeriod as long as step, decompPeriod computes the decomposition and pairlist
decompPeriod ''' + str(steps) + '''
outputFormat dcd

cutoff 34.0
pairlistDistance 500

## Add a dummy particle to make arbd engine happy
particle P
num 1
gridFile ''' + grid_path + '''

## Add rigid body particles
'''
    for count, key in enumerate(diffusible_objects):
        clean_key = key.replace('.aligned', '')
        text += '''
rigidBody ''' + key + '''
num ''' + str(num_copies_per_objects[clean_key]) + '''
mass ''' + str(masses[count]) + '''
inertia ''' + ' '.join([str(elm) for elm in inertias[count]]) + '''
transDamping ''' + ' '.join([str(elm) for elm in transDampings[count]]) + '''
rotDamping ''' + ' '.join([str(elm) for elm in rotDampings[count]]) + '''
densityGrid elec ''' + key + '''.charge.dx
''' + interactive + 'potentiaGrid elec ' + key + '''.elec.smoothed.dx
potentialGridScale elec 0.59616195'''
        for i in range(num_heavy_cluster + 1):
            text += '''
densityGrid vdw''' + str(i) + ' ' + key + '.vdw' + str(i) + '''.den.dx
''' + interactive + 'potentialGrid vdw' + str(i) + ' ' + key + '.vdw' + str(i) + '''.pot.smoothed.dx
potentialGridScale vdw''' + str(i) + ' 0.59616195'
        text += '''

'''

    if len(static_objects) > 0:
        text += '''
# gridFile describes the environement potential. "elec", "vdw1" are labels that match the density of the moving particle to that of the environemnt'''
        for key in static_objects:
            text += '''
gridFile elec ''' + key + '''.elec.smoothed.dx
pmfScale elec 0.59616195'''
            for i in range(num_heavy_cluster + 1):
                text += '''
gridFile vdw''' + str(i) + ' ' + key + '.vdw' + str(i) + '''.pot.smoothed.dx
pmfScale vdw''' + str(i) + ' 0.59616195'

    if len(Extra_pots_tags) > 0:
        for key in Extra_pots_tags:
            dx, vdw = key.replace(' ','').replace('(','').replace(')','').split(',')
            text += '''
gridFile ''' + vdw + ' ' + dx + '''
pmfScale ''' + vdw + ' 0.59616195'

    text += '''
'''

    for key in BDCoordinates:
        text += '''
inputRBCoordinates ''' + key

    with open(out_path, 'w') as fout:
        fout.write(text)

def Write_job_submission_template(ARBD_path='~cmaffeo2/scratch/arbd.dev/src/arbd',
                                 config='Replica_0.bd',
                                 output='Replica_0',
                                 out_path='Replica_0.sh'):
    text = '''#! /bin/bash

''' + ARBD_path + ' -g 0 ' + config + ' ' + output + ' | tee ' + output + '.log'
    with open(out_path, 'w') as fout:
        fout.write(text)

def Write_segmentations(buff=1, out_path='5-segment_Chun.tcl'):
    text = '''#Dissect structure into N pieces
set prefixes $argv
set Inpsf [lindex $prefixes 0]
set Inpdb [lindex $prefixes 1]
set outPrefix [lindex $prefixes 2]
set num_x [lindex $prefixes 3]
set num_y [lindex $prefixes 4]
set num_z [lindex $prefixes 5]
set buff ''' + str(buff) + '''

#Obtain intervals
mol new $Inpsf.psf
mol addfile $Inpdb.pdb

set al [atomselect top "all"]
set mM [measure minmax $al]
set dimAl [vecsub [lindex $mM 1] [lindex $mM 0]]
set min_x [expr [lindex [lindex $mM 0] 0] - $buff]
set min_y [expr [lindex [lindex $mM 0] 1] - $buff]
set min_z [expr [lindex [lindex $mM 0] 2] - $buff]
set dx [expr ([lindex $dimAl 0] + 2 * $buff)/$num_x]
set dy [expr ([lindex $dimAl 1] + 2 * $buff)/$num_y]
set dz [expr ([lindex $dimAl 2] + 2 * $buff)/$num_z]

set count 0
for {set i 0} {$i < $num_x} {incr i} {
  for {set j 0} {$j < $num_y} {incr j} {
      for {set k 0} {$k < $num_z} {incr k} {
        set outName $outPrefix.$count
        set low_x [expr $min_x + $i * $dx]
        set low_y [expr $min_y + $j * $dy]
        set low_z [expr $min_z + $k * $dz]
        set up_x [expr $min_x + ($i+1) * $dx]
        set up_y [expr $min_y + ($j+1) * $dy]
        set up_z [expr $min_z + ($k+1) * $dz]
        set sel [atomselect top "(x > $low_x and x < $up_x) and (y > $low_y and y < $up_y) and (z > $low_z and z < $up_z)"]
        set sel_N [$sel num]
        if {$sel_N > 0} {
          $sel writepqr $outName.pqr
          $sel writepsf $outName.psf
          $sel writepdb $outName.pdb
          set count [expr $count + 1]
        }
    }
  }
}
'''
    with open(out_path, 'w') as fout:
        fout.write(text)

def Write_map_gluing_for_elec(out_path='6-gluing_Chun.tcl'):
    text = '''set prefixes $argv
set InMap1 [lindex $prefixes 0]
set InMap2 [lindex $prefixes 1]
set OutMap [lindex $prefixes 2]

voltool add -i1 $InMap1 -i2 $InMap2 -union -nointerp -o $OutMap
'''
    with open(out_path, 'w') as fout:
        fout.write(text)

def Write_coordinate_files(coors, out_path='test_coor.txt'):
    text = '''# one line per rigid body
# x y z Oxx Oxy Oxz Oyx Oyy Oyz Ozx Ozy Ozz
'''
    for coor in coors:
        text += ' '.join([str(coor[0]), str(coor[1]), str(coor[2])]) + \
        ' 1 0 0 0 1 0 0 0 1\n'

    with open(out_path, 'w') as fout:
        fout.write(text.strip())
