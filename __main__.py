# __main__.py

import sys
import SimpleARBD
from SimpleARBD.SimpleARBD import SetupARBD, SimulateARBD
from SimpleARBD import Read_inputs_for_ARBD
from SimpleARBD import Write_accessory_scripts_for_ARBD

#in_Config = sys.argv[1]

def main():
    in_Config = sys.argv[1]
    if 'in_Config' in locals():
        if SimpleARBD.Read_inputs_for_ARBD.Check_config_name(in_Config):
            Config = Read_inputs_for_ARBD.Read_SimpleARBD_config(in_Config)
            Simple_arbd = SetupARBD(Config)
            print('''..............
Diffusible object(s): {}
.............'''.format(Simple_arbd.diffusible_objects))
            print('''..............
Aligning principle axes of diffusible object(s).
.............''')
            Simple_arbd.Align_diffusible_objs()
            print('''..............
Getting hydrodynamic properties of diffusible object(s).
.............''')
            Simple_arbd.Get_hydrodynamic_properties()
            print('''..............
Getting charge distributions for all object(s).
.............''')
            Simple_arbd.Get_charge_density()
            print('''..............
Getting Eelectrostatic map(s) for all object(s).
.............''')
            Simple_arbd.Get_EM_maps()
            print('''..............
Getting contact force (Vdw) profile(s) for all object(s).
.............''')
            Simple_arbd.Get_VDW_maps()
            print('''..............
Applying whitening to computed represenatations to mimick motion(s) of object(s).
.............''')
            Simple_arbd.Apply_Gaussian_smoothing()
            print('''..............
Reading parameters for ARBD simulations.
.............''')
            Simulate_simple_arbd = SimulateARBD(Simple_arbd, Config)
            print('''..............
Creating boundaries (if any) for ARBD simulations
.............''')
            Simulate_simple_arbd.Create_BD_boundary()
            print('''..............
Generating initial coordinates ARBD simulations.
.............''')
            Simulate_simple_arbd.Generate_initial_coordinates()
            print('''..............
Preparing configuration files for ARBD simulations.
.............''')
            Simulate_simple_arbd.Write_configs()
            print('''..............
Drafting a bash template for running ARBD simulations on a local gpu.
.............''')
            Simulate_simple_arbd.Write_job_submission_template()
            Simulate_simple_arbd.Write_DGX2_submittion_script()
            print('''..............
Transferring files to the destination for ARBD simutlaions.
.............''')
            Simulate_simple_arbd.Tranfer_files()
            print('''..............
Cleaning up the workign directory.
.............''')
            Simulate_simple_arbd.Remove_files()
            pass
        else:
            print('Invalid suffix for the configuration file.')
    else:
        print('No configuration file provided.')

if __name__ == "__main__":
    main()
