#---This is a package for reading files, userinputs, and other miscellaneous for setting up ARBD
import re

def Input_structures_for_diffusible_object():
    structures = []
    quit = False
    print("Please list the path for the structure of each diffusible object involved in your BD simulations (Type Q to quit the input process):")
    count = 1
    while not quit:
        ans_in = input("Path for object {}: ".format(count))
        if ans_in == 'Q':
            quit = True
        else:
            structures.append(ans_in)
    if len(structures) > 0:
        print("The system(s) input are: ", structures)
        #print("Please make sure your objects (psf, pdb) are properly orientated and centered for your study.")
    else:
        print("No object has been input.")
    return structures

def Input_structures_for_static_object():
    structures = []
    quit = False
    print("Please list the path for the structure of each static object involved in your BD simulations (Type Q to quit the input process):")
    count = 1
    while not quit:
        ans_in = input("Path for object {}: ".format(count))
        if ans_in == 'Q':
            quit = True
        else:
            structures.append(ans_in)
    if len(structures) > 0:
        print("The system(s) input are: ", structures)
        print("Please make sure your objects (psf, pdb) are properly orientated and centered for your study.")
    else:
        print("No object has been input.")
    return structures

def Input_salt_concentration():
    ans_in = input("Input the salt concentration (M) for your BD simulation: ")
    try:
        ans_in_num = float(ans_in)
    except:
        ans_in_num = -1
    while ans_in_num < 0:
        print("Enter a positive number.")
        ans_in = input("Input the salt concentration (M) for your BD simulation: ")
        try:
            ans_in_num = float(ans_in)
        except:
            ans_in_num = -1
    return ans_in

def Input_simuation_temperature():
    ans_in = input("Input the temperature (K) for your BD simulation: ")
    try:
        ans_in_num = float(ans_in)
    except:
        ans_in_num = -1
    while ans_in_num < 0:
        print("Enter a positive number.")
        ans_in = input("Input the temperature (K) for your BD simulation: ")
        try:
            ans_in_num = float(ans_in)
        except:
            ans_in_num = -1
    return ans_in

def Read_system_mass(in_path):
    count = 0
    mass = []
    with open(in_path, 'r') as fin:
        for line in fin.readlines():
            in_line = line.strip()
            if len(in_line) > 0:
                mass.append(float(in_line))
                count += 1
            if count == 1:
                break
    return mass[0]

def Read_system_dimension(in_path):
    dimensions = []
    with open(in_path, 'r') as fin:
        for line in fin.readlines():
            in_line = line.strip()
            if len(in_line) > 0:
                dimensions.append(float(in_line))
    return dimensions

def Read_system_inertia(in_path):
    inertias = []
    with open(in_path, 'r') as fin:
        for line in fin.readlines():
            in_line = line.strip()
            if len(in_line) > 0:
                elms = in_line.split(' ')
                if len(elms) == 3:
                    inertias = elms
                    #--should add some data type checking here
                    break
    return inertias

def Read_system_damping_coefficients(in_path):
    tranDampings = []
    rotDampings = []
    count = 0
    with open(in_path, 'r') as fin:
        for line in fin.readlines():
            in_line = line.strip()
            if len(in_line) > 0:
                elms = in_line.split(' ')
                if len(elms) == 3 and count == 0:
                    tranDampings = elms
                    #--should add some data type checking here
                    count += 1
                    continue
                if len(elms) == 3 and count == 1:
                    rotDampings = elms
                    #--should add some data type checking here
                    count += 1
                    break
    return (tranDampings, rotDampings)

def Read_SimpleARBD_config(in_path):
    Config = {}
    with open(in_path, 'r') as fin:
        text = fin.read()
    m = (re.search(r'Diffusible_objects:([ \w]+)', text))
    Config['diffusible_objects'] = m.group(1).strip().split(' ')
    m = (re.search(r'Static_objects \(Enter NA for no static object\):([ \w]+)', text))
    if m.group(1).strip() == 'NA':
        Config['static_objects'] = []
    else:
        Config['static_objects'] = m.group(1).strip().split(' ')
    m = (re.search(r'SaltConcentration:(\s*[0-9]*.[0-9]*)', text))
    Config['saltConcentration'] = m.group(1).strip()
    m = (re.search(r'Temperature \(K\):(\s*[0-9]*.*[0-9]*)', text))
    Config['temperature'] = m.group(1).strip()
    m = (re.search(r'Viscosity:(\s*[0-9]*.*[0-9]*)', text))
    Config['viscosity'] = m.group(1).strip()
    m = (re.search(r'Solvent_density:(\s*[0-9]*.*[0-9]*)', text))
    Config['solvent_density'] = m.group(1).strip()
    m = (re.search(r'Number_of_heavy_cluster \(Integer\):(\s*[0-9]+)', text))
    Config['num_heavy_cluster'] = m.group(1).strip()
    m = (re.search(r'GaussianWidth:(\s*[0-9]*.*[0-9]*)', text))
    Config['GaussianWidth'] = m.group(1).strip()
    m = (re.search(r'Skip_parametrizing_diffusible \(Yes/No\):([ \w]+)', text))
    Config['Skip_parametrizing_diffusible'] = m.group(1).strip()
    m = (re.search(r'Gigantic_stat_objects \(Yes/No\):([ \w]+)', text))
    Config['Gigantic_stat_objects'] = m.group(1).strip()
    m = (re.search(r'Python_path:(\s*\S+)', text))
    Config['python_path'] = m.group(1).strip()
    m = (re.search(r'Hydro_path:(\s*\S+)', text))
    Config['hydro_path'] = m.group(1).strip()
    m = (re.search(r'Apbs_path:(\s*\S+)', text))
    Config['apbs_path'] = m.group(1).strip()
    m = (re.search(r'Vmd_path:(\s*\S+)', text))
    Config['vmd_path'] = m.group(1).strip()
    m = (re.search(r'Parameters_folder:(\s*\S+)', text))
    Config['parameters_folder'] = m.group(1).strip()
    m = (re.search(r'Num_replicas \(Integer\):(\s*[0-9]+)', text))
    Config['num_replicas'] = m.group(1).strip()
    m = (re.search(r'Timestep \(Float\):(\s*[0-9]*.*[0-9]*)', text))
    Config['timestep'] = m.group(1).strip()
    m = (re.search(r'Steps \(Integer\):(\s*[0-9]+)', text))
    Config['steps'] = m.group(1).strip()
    m = (re.search(r'Interactive \(Yes/No\):([ \w]+)', text))
    if m.group(1) == 'Yes':
        Config['interactive'] = ''
    else:
        Config['interactive'] = '#'
    m = (re.search(r'Grid_path:(\s*\S+)', text))
    Config['grid_path'] = m.group(1).strip()
    m = (re.search(r'Number_of_copies_per_object \(Integer\(s\)\):([ 0-9]+)', text))
    Config['num_copies_per_object'] = m.group(1).strip().split(' ')

    temp = re.search(r'(Extra_potentials_tags \(Path, vdw cluster group\):[\s\S]*)\n', text)
    m = re.findall(r'\((\S+\.dx,\s*\w+)\)', temp.group(0))
    Config['Extra_pots_tags'] = m

    m = (re.search(r'WellDepth \(Positive\):\s*([0-9]+[\.]*[0-9]*)', text))
    Config['wellDepth'] = m.group(1).strip()
    m = (re.search(r'WellResolution \(Positive\):\s*([0-9]+[\.]*[0-9]*)', text))
    Config['wellResolution'] = m.group(1).strip()
    m = (re.search(r'CellBasisVector1:\s*([0-9]+[\.]*[0-9]* [0-9]+[\.]*[0-9]* [0-9]+[\.]*[0-9]*)', text))
    Config['cellBasisVector1'] = m.group(1).strip()
    m = (re.search(r'CellBasisVector2:\s*([0-9]+[\.]*[0-9]* [0-9]+[\.]*[0-9]* [0-9]+[\.]*[0-9]*)', text))
    Config['cellBasisVector2'] = m.group(1).strip()
    m = (re.search(r'CellBasisVector3:\s*([0-9]+[\.]*[0-9]* [0-9]+[\.]*[0-9]* [0-9]+[\.]*[0-9]*)', text))
    Config['cellBasisVector3'] = m.group(1).strip()
    m = (re.search(r'CellOrigin:\s*(-*[0-9]+[\.]*[0-9]* -*[0-9]+[\.]*[0-9]* -*[0-9]+[\.]*[0-9]*)', text))
    Config['cellOrigin'] = m.group(1).strip()
    m = (re.search(r'InitialCoorBasisVector1:\s*([0-9]+[\.]*[0-9]* [0-9]+[\.]*[0-9]* [0-9]+[\.]*[0-9]*)', text))
    Config['InitialCoorBasisVector1'] = m.group(1).strip()
    m = (re.search(r'InitialCoorBasisVector2:\s*([0-9]+[\.]*[0-9]* [0-9]+[\.]*[0-9]* [0-9]+[\.]*[0-9]*)', text))
    Config['InitialCoorBasisVector2'] = m.group(1).strip()
    m = (re.search(r'InitialCoorBasisVector3:\s*([0-9]+[\.]*[0-9]* [0-9]+[\.]*[0-9]* [0-9]+[\.]*[0-9]*)', text))
    Config['InitialCoorBasisVector3'] = m.group(1).strip()
    m = (re.search(r'InitialCoorOrigin:\s*(-*[0-9]+[\.]*[0-9]* -*[0-9]+[\.]*[0-9]* -*[0-9]+[\.]*[0-9]*)', text))
    Config['InitialCoorOrigin'] = m.group(1).strip()
    #m = (re.search(r'RBCoordinates:([ \S]+)', text))
    #Config['RBCoordinates'] = m.group(1).strip().split(' ')
    m = (re.search(r'ARBD_path:(\s*\S+)', text))
    Config['ARBD_path'] = m.group(1).strip()

    #---read the path for simulation
    m = re.search(r'Path_for_ARBD_simulations:(\s*\S+)', text)
    Config['Path_for_ARBD_simulations'] = m.group(1).strip()

    return Config

def Check_config_name(in_Config):
    m = re.search(r'(\S+).SimpleARBD_conf$', in_Config)
    if m:
        print('Parsing: ', in_Config)
        return True
    return False
