# /software/python3-3.4.1/bin/python3

import os
import shutil
import glob
from os.path import isfile
from SimpleARBD import Read_inputs_for_ARBD
from SimpleARBD import Write_accessory_scripts_for_ARBD
from SimpleARBD import Accessory_routines_for_ARBD

class SetupARBD(object):

    def __init__(self, Config):

        self.diffusible_objects = Config['diffusible_objects']
        self.static_objects = Config['static_objects']
        if len(self.static_objects) > 0:
            systems = [key + '.aligned' for key in self.diffusible_objects] + \
            self.static_objects
        else:
            systems = [key + '.aligned' for key in self.diffusible_objects]
        self.computed_systems = systems
        self.saltConcentration = float(Config['saltConcentration'])
        self.temperature = float(Config['temperature'])
        self.viscosity = float(Config['viscosity'])
        self.solvent_density = float(Config['solvent_density'])
        self.number_of_heavy_atom_clusters = int(Config['num_heavy_cluster'])
        self.GaussianWidth = float(Config['GaussianWidth'])
        self.Skip_parametrizing_diffusible = str(Config['Skip_parametrizing_diffusible'])
        self.Gigantic_stat_objects = str(Config['Gigantic_stat_objects'])
        self.python_path = Config['python_path']
        self.hydro_path = Config['hydro_path']
        self.apbs_path = Config['apbs_path']
        self.vmd_path = Config['vmd_path']
        self.parameters_folder = Config['parameters_folder']

    def Align_diffusible_objs(self, tcl_path='1-align_Chun.tcl'):
        Write_accessory_scripts_for_ARBD.Write_align_tcl(out_path=tcl_path)
        if self.Skip_parametrizing_diffusible in ["No", "no"]:
            for key in self.diffusible_objects:
                cmd_in = self.vmd_path + ' -dispdev text -args ' + key + ' < ' + tcl_path
                os.system(cmd_in)

    def Get_hydrodynamic_properties(self):
        if self.Skip_parametrizing_diffusible in ["No", "no"]:
            temperature = self.temperature - 273
            for key in self.diffusible_objects:
                mass_path = key + '.mass.txt'
                mass = Read_inputs_for_ARBD.Read_system_mass(mass_path)
                system = key + '.aligned'
                Write_accessory_scripts_for_ARBD.Write_hydropro_config(
                    key, mass, temperature=temperature,
                    viscosity=self.viscosity,
                    solvent_density=self.solvent_density
                )
                cmd_in = '''(
                        cd ./
                        ln -b -s ''' + system + '''.pdb hydro.pdb
                        ''' + self.hydro_path + '''
                        )'''
                os.system(cmd_in)
                os.system('unlink hydro.pdb')
                hydroproFile = key + '.hydro-res.txt'
                massFile = key + '.mass.txt'
                inertiaFile = key + '.inertia.txt'
                outFile =  key + '.damping-coeffs.txt'
                Accessory_routines_for_ARBD.Get_damping_coefficients(
                    hydroproFile, massFile, inertiaFile, outFile
                )


    def Get_charge_density(self,resolution=2, tcl_path="2-charge-density_Chun.tcl"):
        Write_accessory_scripts_for_ARBD.Write_charge_density_tcl(
            resolution=resolution,
            out_path=tcl_path
        )
        if self.Skip_parametrizing_diffusible in ["No", "no"]:
            for key in self.computed_systems:
                cmd_in = self.vmd_path + ' -dispdev text -args ' + key + ' < ' + tcl_path
                os.system(cmd_in)
                inFile = key + '.chargeDensity.dx'

                os.system('rm -f temp0.dx')
                os.system('rm -f temp1.dx')

                outFile = key + '.charge.dx'
                netChargeFile = key + '.netCharge.dat'
                Accessory_routines_for_ARBD.Fix_charge(inFile, outFile, netChargeFile)
        elif self.Skip_parametrizing_diffusible in ["Yes", "yes"]:
            for key in self.static_objects:
                cmd_in = self.vmd_path + ' -dispdev text -args ' + key + ' < ' + tcl_path
                os.system(cmd_in)
                inFile = key + '.chargeDensity.dx'
                outFile = key + '.charge.dx'
                netChargeFile = key + '.netCharge.dat'
                Accessory_routines_for_ARBD.Fix_charge(inFile, outFile, netChargeFile)


    def Get_EM_maps(self, buffer=50):
        systems = [key + '.aligned' for key in self.diffusible_objects]
        if self.Skip_parametrizing_diffusible in ["No", "no"]:
            for key in systems:
                in_path = key + '.dimension.dat'
                dimensions = Read_inputs_for_ARBD.Read_system_dimension(in_path)
                Write_accessory_scripts_for_ARBD.Write_apbs_config(
                    system=key, in_xyz_dim=dimensions, conc=self.saltConcentration,
                    temperature=self.temperature, buffer=buffer
                )
                cmd_in = self.apbs_path + ' ' + key + '.apbs'
                os.system(cmd_in)
                inFile = key + '.elec.tmp.dx'
                outFile = key + '.elec.dx'
                Accessory_routines_for_ARBD.Bound_grid(inFile, outFile, -20, 20)

        if self.Gigantic_stat_objects in ["No", "no"]:
            Large_system = 'Off'
            if len(self.static_objects) > 0:
                for key in self.static_objects:
                    in_path = key + '.dimension.dat'
                    dimensions = Read_inputs_for_ARBD.Read_system_dimension(in_path)
                    Write_accessory_scripts_for_ARBD.Write_apbs_config(
                        system=key, in_xyz_dim=dimensions, conc=self.saltConcentration,
                        temperature=self.temperature, buffer=buffer, Large_system=Large_system
                    )
                    cmd_in = self.apbs_path + ' ' + key + '.apbs'
                    os.system(cmd_in)
                    inFile = key + '.elec.tmp.dx'
                    outFile = key + '.elec.dx'
                    Accessory_routines_for_ARBD.Bound_grid(inFile, outFile, -20, 20)
        elif self.Gigantic_stat_objects in ["Yes", "yes"]:
            Large_system = 'On'
            if len(self.static_objects) > 0:
                self.Gigantic_stat_objects_segmentation_num = {}
                for key in self.static_objects:
                    in_path = key + '.dimension.dat'
                    tcl_path = '5-segment_Chun.tcl'
                    out_temp = key + '.stat_temp'
                    dimensions = Read_inputs_for_ARBD.Read_system_dimension(in_path)
                    threshold=300
                    nx, ny, nz = Accessory_routines_for_ARBD.Find_segments_num(dimensions, threshold=threshold)
                    Write_accessory_scripts_for_ARBD.Write_segmentations(out_path=tcl_path)
                    cmd_in = ' '.join([
                        'VMDNOCUDA=1', self.vmd_path, '-dispdev text -args', key, key, out_temp, str(nx), str(ny), str(nz), '<', tcl_path
                        ])
                    os.system(cmd_in)
                    tcl_path = '6-gluing_Chun.tcl'
                    Write_accessory_scripts_for_ARBD.Write_map_gluing_for_elec(out_path=tcl_path)
                    count = 0
                    for i in range(nx * ny * nz):
                        temp_sys = out_temp + '.' + str(i)
                        pqr_file = temp_sys + '.pqr'
                        if isfile(pqr_file):
                            dimensions = [threshold] * 3
                            Write_accessory_scripts_for_ARBD.Write_apbs_config(
                                system=temp_sys, in_xyz_dim=dimensions, conc=self.saltConcentration,
                                temperature=self.temperature, buffer=buffer, Large_system=Large_system
                            )
                            cmd_in = self.apbs_path + ' ' + temp_sys + '.apbs'
                            os.system(cmd_in)
                            if i == 0:
                                inFile = temp_sys + '.elec.tmp.dx'
                                outFile = temp_sys + '.elec.dx'
                                os.system(' '.join(['cp', inFile, outFile]))
                            elif i > 0:
                                count += 1
                                inFile1 = out_temp + '.' + str(i-1) + '.elec.dx'
                                inFile2 = temp_sys + '.elec.tmp.dx'
                                outFile = temp_sys + '.elec.dx'
                                cmd_in = ' '.join([
                                    'VMDNOCUDA=1', self.vmd_path, '-dispdev text -args', inFile1, inFile2, outFile, '<', tcl_path
                                ])
                                os.system(cmd_in)
                    self.Gigantic_stat_objects_segmentation_num[key] = count

                    inFile = out_temp + '.' + str(count) + '.elec.dx'
                    if isfile(inFile):
                        outFile = key + '.elec.tmp.dx'
                        os.system(' '.join(['cp', inFile, outFile]))
                        inFile = outFile
                        outFile = key + '.elec.dx'
                        Accessory_routines_for_ARBD.Bound_grid(inFile, outFile, -20, 20)


    def Get_VDW_maps(self, potResolution=1, denResolution=2, tcl_path="3-vdw_Chun.tcl"):

        num_cluster = self.number_of_heavy_atom_clusters
        if self.Skip_parametrizing_diffusible in ["No", "no"]:
            for key in self.diffusible_objects:
                system = key + '.aligned'
                Write_accessory_scripts_for_ARBD.Write_vdw_den_pot_tcl(
                    potResolution=potResolution,
                    denResolution=denResolution,
                    num_heavy_cluster=self.number_of_heavy_atom_clusters,
                    python_path=self.python_path,
                    out_path=tcl_path,
                    parameter_folder=self.parameters_folder,
                    diffuse_or_static='diffuse'
                )
                cmd_in = 'VMDNOCUDA=1 ' + self.vmd_path + ' -dispdev text -args ' + system + ' < ' + tcl_path
                os.system(cmd_in)


        if len(self.static_objects) > 0:
            if self.Gigantic_stat_objects in ["No", "no"]:
                for key in self.static_objects:
                    Write_accessory_scripts_for_ARBD.Write_vdw_den_pot_tcl(
                        potResolution=potResolution,
                        denResolution=denResolution,
                        num_heavy_cluster=self.number_of_heavy_atom_clusters,
                        python_path=self.python_path,
                        out_path=tcl_path,
                        parameter_folder=self.parameters_folder,
                        diffuse_or_static='static'
                    )
                    cmd_in = 'VMDNOCUDA=1 ' + self.vmd_path + ' -dispdev text -args ' + key + ' < ' + tcl_path
                    os.system(cmd_in)

            elif self.Gigantic_stat_objects in ["Yes", "yes"]:
                Write_accessory_scripts_for_ARBD.Write_vdw_den_pot_tcl(
                    potResolution=2*potResolution,
                    denResolution=denResolution,
                    num_heavy_cluster=self.number_of_heavy_atom_clusters,
                    python_path=self.python_path,
                    out_path=tcl_path,
                    parameter_folder=self.parameters_folder,
                    diffuse_or_static='static'
                )
                for key in self.static_objects:
                    segment_num = self.Gigantic_stat_objects_segmentation_num[key]
                    systems_temp = [key + '.stat_temp.' + str(count) for count in range(segment_num + 1)]
                    for system in systems_temp:
                        cmd_in = 'VMDNOCUDA=1 ' + self.vmd_path + ' -dispdev text -args ' + system + ' < ' + tcl_path
                        os.system(cmd_in)

                    #--Gluing maps
                    tcl_path = '6-gluing_Chun.tcl'
                    for j in range(segment_num + 1):
                        for i in range(num_cluster + 1):
                            if j == 0:
                                inFile = key + '.stat_temp.' + str(j) + '.vdw' + str(i) + '.pot.dx'
                                outFile = key + '.stat_temp.' + str(j) + '.vdw' + str(i) + '.pot.tmp.dx'
                                os.system(' '.join(['cp', inFile, outFile]))
                            elif j > 0:
                                inFile1 = key + '.stat_temp.' + str(j-1) + '.vdw' + str(i) + '.pot.tmp.dx'
                                inFile2 = key + '.stat_temp.' + str(j) + '.vdw' + str(i) + '.pot.dx'
                                outFile = key + '.stat_temp.' + str(j) + '.vdw' + str(i) + '.pot.tmp.dx'
                                cmd_in = ' '.join([
                                    'VMDNOCUDA=1', self.vmd_path, '-dispdev text -args', inFile1, inFile2, outFile, '<', tcl_path
                                ])
                                os.system(cmd_in)

                    for i in range(num_cluster + 1):
                        inFile = key + '.stat_temp.' + str(j) + '.vdw' + str(i) + '.pot.tmp.dx'
                        if isfile(inFile):
                            outFile = key + '.vdw' + str(i) + '.pot.tmp.dx'
                            os.system(' '.join(['cp', inFile, outFile]))
                            inFile = outFile
                            outFile = key + '.vdw' + str(i) + '.pot.dx'
                            Accessory_routines_for_ARBD.Bound_grid(inFile, outFile, -20, 20)



    def Apply_Gaussian_smoothing(self, tcl_path='4-smooth_Chun.tcl'):
        for key in self.computed_systems:
            #--EM
            in_file = key + '.elec.dx'
            out_file = key + '.elec.smoothed.dx'
            Write_accessory_scripts_for_ARBD.Write_smoothing_tcl(
                in_file=in_file,
                out_file=out_file,
                GaussianWidth=self.GaussianWidth,
                out_path=tcl_path
            )
            cmd_in = 'vmd -dispdev text < ' + tcl_path
            os.system(cmd_in)

            #--Vdw
            for i in range(self.number_of_heavy_atom_clusters + 1):
                in_file = key + '.vdw' + str(i) + '.pot.dx'
                out_file = key + '.vdw' + str(i) + '.pot.smoothed.dx'
                Write_accessory_scripts_for_ARBD.Write_smoothing_tcl(
                    in_file=in_file,
                    out_file=out_file,
                    GaussianWidth=self.GaussianWidth,
                    out_path=tcl_path
                )
                cmd_in = 'vmd -dispdev text < ' + tcl_path
                os.system(cmd_in)

class SimulateARBD(object):

    def __init__(self, ARBD_object, Config):
        self.ARBD_object = ARBD_object
        self.num_replicas = int(Config['num_replicas'])
        self.timestep = float(Config['timestep'])
        self.steps = int(Config['steps'])
        self.interactive = Config['interactive']
        self.grid_path = Config['grid_path']
        self.num_copies_per_object = {
        elm[1]: int(Config['num_copies_per_object'][elm[0]]) for elm in enumerate(self.ARBD_object.diffusible_objects)
        }
        self.Extra_pots_tags = Config['Extra_pots_tags']
        self.wellDepth = -1 * float(Config['wellDepth'])
        self.wellResolution = float(Config['wellResolution'])
        if len(Config['Extra_pots_tags']) > 0:
            self.use_Boundary = True
            self.Boundary_path = Config['Extra_pots_tags'][0].replace(',', '').split()[0]
            self.Boundary_vdwType = Config['Extra_pots_tags'][0].replace(',', '').split()[1]
            self.cellBasisVector1 = [round(float(elm)) for elm in Config['cellBasisVector1'].split()]
            self.cellBasisVector2 = [round(float(elm)) for elm in Config['cellBasisVector2'].split()]
            self.cellBasisVector3 = [round(float(elm)) for elm in Config['cellBasisVector3'].split()]
            self.cellOrigin = [round(float(elm)) for elm in Config['cellOrigin'].split()]
        #self.RBCoordinates = Config['RBCoordinates']
        self.InitialCoorBasisVector1 = [round(float(elm)) for elm in Config['InitialCoorBasisVector1'].split()]
        self.InitialCoorBasisVector2 = [round(float(elm)) for elm in Config['InitialCoorBasisVector2'].split()]
        self.InitialCoorBasisVector3 = [round(float(elm)) for elm in Config['InitialCoorBasisVector3'].split()]
        self.InitialCoorOrigin = [round(float(elm)) for elm in Config['InitialCoorOrigin'].split()]
        self.ARBD_path = Config['ARBD_path']
        self.Simulation_path = Config['Path_for_ARBD_simulations']

    def Create_BD_boundary(self, resolution=2, blur=5, wellDepth = -1, scale = 1.5):

        if self.use_Boundary:
            wellDepth = self.wellDepth
            resolution = self.wellResolution
            v1 = self.cellBasisVector1
            v2 = self.cellBasisVector2
            v3 = self.cellBasisVector3
            out_path = self.Boundary_path
            cellOrigin = self.cellOrigin
            #--search for possible resolutions
            dd, n1, n2, n3 = Accessory_routines_for_ARBD.Find_boundary_resolution(v1, v2, v3, resolution)

            #--Create a rectangular mesh
            mx, my, mz, Mx, My, Mz = Accessory_routines_for_ARBD.Find_boundary_end_points(v1, v2, v3, cellOrigin)
            wallRangeX = (scale*(mx+2*blur), scale*(Mx+2*blur))
            wallRangeY = (scale*(my+2*blur), scale*(My+2*blur))
            wallRangeZ = (scale*(mz+2*blur), scale*(Mz+2*blur))
            X, Y, Z, dx, dy, dz = Accessory_routines_for_ARBD.Create_a_rectangular_mesh(
            wallRangeX, wallRangeY, wallRangeZ, blur, dd
            )

            #--Define the Boundary region
            region_pts = Accessory_routines_for_ARBD.Create_boundary_region(
                mx, my, mz, n1, n2, n3, cellBasisVector1=v1, cellBasisVector2=v2, cellBasisVector3=v3
            )
            mesh_pts = Accessory_routines_for_ARBD.Convert_math_pt_to_mesh_pt(
                region_pts, wallRangeX[0], wallRangeY[0], wallRangeZ[0], dx, dy, dz
            )

            #--Create the boundary
            Accessory_routines_for_ARBD.Create_the_well(
            X, wellDepth, mesh_pts, blur, dx, dy, dz, wallRangeX[0], wallRangeY[0], wallRangeZ[0], out_path
            )

    def Generate_initial_coordinates(self):

        for key in self.ARBD_object.diffusible_objects:
            #--Step 1: Read the dimension of a system
            in_path = key + '.aligned.dimension.dat'
            dimensions = Read_inputs_for_ARBD.Read_system_dimension(in_path)

            #--Step 2: Generate the basis vectors for spanning the initial coordinates for the diffusibe objects
            bv1, bv2, bv3, n1, n2, n3 = Accessory_routines_for_ARBD.Generate_spanning_vectors(
            self.InitialCoorBasisVector1, self.InitialCoorBasisVector2, self.InitialCoorBasisVector3, dimensions
            )

            #--Step 3: Generate a list of coordinates for each replica
            for replica_index in range(self.num_replicas):
                coors = Accessory_routines_for_ARBD.Generate_coordinates(
                bv1, bv2, bv3, n1, n2, n3, self.num_copies_per_object[key], self.InitialCoorOrigin, replica_index
                )
                #--STep 4: Write out coordinates
                out_path = key + '_bdcoor' + str(replica_index) + '.txt'
                Write_accessory_scripts_for_ARBD.Write_coordinate_files(coors, out_path)


    def Write_configs(self):

        grid_path = self.grid_path
        if not isfile(grid_path):
            Accessory_routines_for_ARBD.Create_null(grid_path=grid_path)

        masses = [
           Read_inputs_for_ARBD.Read_system_mass(key + '.mass.txt') for key in self.ARBD_object.diffusible_objects
        ]
        inertias = [
           Read_inputs_for_ARBD.Read_system_inertia(key + '.inertia.txt') for key in self.ARBD_object.diffusible_objects
        ]
        Dampings = [
           Read_inputs_for_ARBD.Read_system_damping_coefficients(key + '.damping-coeffs.txt') for key in self.ARBD_object.diffusible_objects
        ]
        transDampings = [Damping[0] for Damping in Dampings]
        rotDampings = [Damping[1] for Damping in Dampings]


        for i in range(self.num_replicas):
            out_path = 'Replica_' + str(i) + '.bd'
            coordinates = [str(key) + '_bdcoor' + str(i) + '.txt' for key in self.ARBD_object.diffusible_objects]
            Write_accessory_scripts_for_ARBD.Write_BD_configs(
                timestep=self.timestep,
                steps=self.steps,
                temperature=self.ARBD_object.temperature,
                num_heavy_cluster=self.ARBD_object.number_of_heavy_atom_clusters,
                interactive=self.interactive,
                grid_path=self.grid_path,
                diffusible_objects=[key + '.aligned' for key in self.ARBD_object.diffusible_objects],
                num_copies_per_objects=self.num_copies_per_object,
                masses=masses,
                inertias=inertias,
                transDampings=transDampings,
                rotDampings=rotDampings,
                static_objects=[key for key in self.ARBD_object.static_objects],
                Extra_pots_tags=self.Extra_pots_tags,
                BDCoordinates=coordinates,
                out_path=out_path
            )

    def Write_job_submission_template(self):
        Write_accessory_scripts_for_ARBD.Write_job_submission_template(
            ARBD_path=self.ARBD_path)

    def Tranfer_files(self):
        current_path = os.getcwd()
        #--check if the diretcory exists
        if os.path.exists(self.Simulation_path):
            print("{} exists. Moving files from {} to {} now.".format(self.Simulation_path, current_path, self.Simulation_path))
        else:
            print("{} does not exist. Creating directory and moving files from {} to {} now.".format(self.Simulation_path, current_path, self.Simulation_path))
            os.makedirs(self.Simulation_path)

        #--move extra potential
        count = 0
        if len(self.Extra_pots_tags) > 0:
            for extra_pot in self.Extra_pots_tags:
                in_path = self.Extra_pots_tags[0].replace(',','').split()[0]
                in_file = os.path.basename(in_path)
                out_path = os.path.join(self.Simulation_path, in_file)
                if count == 0:
                    in_path = os.path.join(current_path, in_file)
                if isfile(in_path):
                    shutil.move(in_path, out_path)
                    print('Moved: ', in_file)
                count += 1

        #--transfer null.dx
        in_path = self.grid_path
        in_file = os.path.basename(in_path)
        in_dir = os.path.dirname(in_path)
        out_path = os.path.join(self.Simulation_path, in_file)
        if not in_dir:
            in_path = os.path.join(current_path, in_file)
        if isfile(in_path):
            shutil.move(in_path, out_path)
            print('Moved: ', in_file)

        #--transfer smoothed dx
        pattern = "/*smoothed*dx"
        in_paths = glob.glob(current_path + pattern)
        for in_path in in_paths:
            in_file = os.path.basename(in_path)
            out_path = os.path.join(self.Simulation_path, in_file)
            if isfile(in_path):
                shutil.move(in_path, out_path)
                print('Moved: ', in_file)

        #--transfer density files
        pattern = "/*den*dx"
        in_paths = glob.glob(current_path + pattern)
        for in_path in in_paths:
            in_file = os.path.basename(in_path)
            out_path = os.path.join(self.Simulation_path, in_file)
            if isfile(in_path):
                shutil.move(in_path, out_path)
                print('Moved: ', in_file)

        #--transfer charge files
        pattern = "/*charge*dx"
        in_paths = glob.glob(current_path + pattern)
        for in_path in in_paths:
            in_file = os.path.basename(in_path)
            out_path = os.path.join(self.Simulation_path, in_file)
            if isfile(in_path):
                shutil.move(in_path, out_path)
                print('Moved: ', in_file)

        #--transfer coordinate files & other useful txt fies
        pattern = "/*txt"
        in_paths = glob.glob(current_path + pattern)
        for in_path in in_paths:
            in_file = os.path.basename(in_path)
            out_path = os.path.join(self.Simulation_path, in_file)
            if isfile(in_path):
                shutil.move(in_path, out_path)
                print('Moved: ', in_file)

        #--transfer configuration files
        pattern = "/*bd"
        in_paths = glob.glob(current_path + pattern)
        for in_path in in_paths:
            in_file = os.path.basename(in_path)
            out_path = os.path.join(self.Simulation_path, in_file)
            if isfile(in_path):
                shutil.move(in_path, out_path)
                print('Moved: ', in_file)

        #--transfer useful dx files
        for key in self.ARBD_object.computed_systems:
            pattern = "/" + key + '.vdw*.pot.dx'
            in_paths = glob.glob(current_path + pattern)
            for in_path in in_paths:
                in_file = os.path.basename(in_path)
                out_path = os.path.join(self.Simulation_path, in_file)
                if isfile(in_path):
                    shutil.move(in_path, out_path)
                    print('Moved: ', in_file)

            in_path = key + '.elec.dx'
            in_file = os.path.basename(in_path)
            out_path = os.path.join(self.Simulation_path, in_file)
            in_dir = os.path.dirname(in_path)
            if not in_dir:
                in_path = os.path.join(current_path, in_file)
            if isfile(in_path):
                shutil.move(in_path, out_path)
                print('Moved: ', in_file)

        #--transfer initial structures
        for key in self.ARBD_object.computed_systems:
            for suffix in [".psf", ".pdb",".pqr"]:
                in_path = key + suffix
                in_file = os.path.basename(in_path)
                out_path = os.path.join(self.Simulation_path, in_file)
                in_dir = os.path.dirname(in_path)
                if not in_dir:
                    in_path = os.path.join(current_path, in_file)
                if isfile(in_path):
                    shutil.move(in_path, out_path)
                    print('Moved: ', in_file)


        #--move dat files
        pattern = "/*dat"
        in_paths = glob.glob(current_path + pattern)
        for in_path in in_paths:
            in_file = os.path.basename(in_path)
            out_path = os.path.join(self.Simulation_path, in_file)
            if isfile(in_path):
                shutil.move(in_path, out_path)
                print('Moved: ', in_file)

        #--move hydro files
        pattern = "/*.hydro*"
        in_paths = glob.glob(current_path + pattern)
        for in_path in in_paths:
            in_file = os.path.basename(in_path)
            out_path = os.path.join(self.Simulation_path, in_file)
            if isfile(in_path):
                shutil.move(in_path, out_path)
                print('Moved: ', in_file)

        #--move submission script
        #--move dat files
        pattern = "/*sh"
        in_paths = glob.glob(current_path + pattern)
        for in_path in in_paths:
            in_file = os.path.basename(in_path)
            out_path = os.path.join(self.Simulation_path, in_file)
            if isfile(in_path):
                shutil.move(in_path, out_path)
                print('Moved: ', in_file)

    def Remove_files(self):
        current_path = os.getcwd()

        pattern = "/*temp*dx"
        in_paths = glob.glob(current_path + pattern)
        for in_path in in_paths:
            os.remove(in_path)

        pattern = "/*stat_temp*p*"
        in_paths = glob.glob(current_path + pattern)
        for in_path in in_paths:
            os.remove(in_path)

        pattern = "/*tcl"
        in_paths = glob.glob(current_path + pattern)
        for in_path in in_paths:
            os.remove(in_path)

        pattern = "/*tmp.dx"
        in_paths = glob.glob(current_path + pattern)
        for in_path in in_paths:
            os.remove(in_path)

        #pattern = "/fort*"
        #in_paths = glob.glob(current_path + pattern)
        #for in_path in in_paths:
        #    os.remove(in_path)

        #pattern = "/io.mc"
        #in_paths = glob.glob(current_path + pattern)
        #for in_path in in_paths:
        #    os.remove(in_path)

        pattern = "/*apbs"
        in_paths = glob.glob(current_path + pattern)
        for in_path in in_paths:
            os.remove(in_path)

    def Write_DGX2_submittion_script(self, out_path='DGX2_script.sh'):
        text = '''#! /bin/bash

#qsub -q dgx -S /bin/bash
#qsub -q dgx -hold_jid ${JOBID} -S /bin/bash

cd ''' + self.Simulation_path + '''

ARBD=~cmaffeo2/scratch/arbd.dev/src/arbd

i=0
for f in *.bd; do
    pre=${f%.bd}
    $ARBD -g $i $f $pre | tee $pre.log &
    sleep 1.5s
    i=$(($i+1))
done

wait'''

        with open(out_path, 'w') as fout:
            fout.write(text)
