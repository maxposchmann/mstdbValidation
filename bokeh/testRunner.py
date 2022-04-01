import json
import parseTests
import subprocess

def runVaporPressures(series):
    inputScript = 'vaporPressure.ti'
    if series['database'] == 'fluoride':
        datapath = fluoridepath
    elif series['database'] == 'chloride':
        datapath = chloridepath
    else:
        print('database not recognized')
        return
    elements = list(series['composition'].keys())
    nElements = len(elements)
    compString = " ".join([str(series['composition'][element]) for element in elements])
    logErr = series['logarithmic error']
    # Write input script
    with open(inputScript, 'w') as inputFile:
        inputFile.write('! Python-generated input file for Thermochimica\n')
        inputFile.write(f'data file         = {datapath}\n')
        inputFile.write(f'temperature unit  = {tunit}\n')
        inputFile.write(f'pressure unit     = {punit}\n')
        inputFile.write(f'mass unit         = {munit}\n')
        inputFile.write(f'nEl               = {nElements}\n')
        inputFile.write(f'iEl               = {" ".join([str(atomic_number_map.index(element)+1) for element in elements])}\n')
        inputFile.write(f'nCalc             = {len(list(currentSeries["samples"].keys()))}\n')
        for sample in currentSeries['samples']:
            s = currentSeries['samples'][sample]
            inputFile.write(f'{s["temperature"]} {press} {compString}\n')
    subprocess.run([thermopath,inputScript])

# Set file names for input/output
infilename  = 'verificationData.json'
outfilename = 'verificationData-tested.json'

# Set path to Thermochimica and outputs
thermopath  = '/media/max/data/thermochimicastuff/thermochimica/bin/RunCalculationList'
outjsonpath = '/media/max/data/thermochimicastuff/thermochimica/thermoout.json'

# Set MSTDB database paths
fluoridepath = '/media/max/data/mstdbValidation/mstdb/Models and Documentation/MSTDB-TC_V1.2_Fluorides_8-0.dat'
chloridepath = '/media/max/data/mstdbValidation/mstdb/Models and Documentation/MSTDB-TC_V1.2_Chlorides_8-0.dat'

# Set units
tunit = 'K'
punit = 'atm'
munit = 'moles'
# Assume all pressures will be 1 atm for now
press = 1.0

atomic_number_map = [
    'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P',
    'S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
    'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh',
    'Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd',
    'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re',
    'Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
    'Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db',
    'Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts', 'Og'
]

# Get data
data = parseTests.jsonTestData(infilename)
nExpRef = len(data.experimentalReferences)
nTestSeries = len(data.testSeries)

for sourceName in data.data['sources']:
    source = data.data['sources'][sourceName]
    # Loop over all associated test series and execture
    for testType in source['tests']:
        for series in source['tests'][testType]:
            print(series)
            currentSeries = source['tests'][testType][series]
            if 'enabled' in currentSeries.keys():
                if currentSeries['enabled'] in ['false', 'False', 'FALSE', 0]:
                    print('disabled')
                    continue
            if testType in ['vapor pressures']:
                runVaporPressures(currentSeries)
                for sample in currentSeries['samples']:
                    print(f"{sample}: {currentSeries['samples'][sample]['status']}")
            # elif testType == 'phase transitions':
            #     runPhaseTransitions(currentSeries)
            # elif testType in ['solubility limits']:
            #     runSolubilityLimits(currentSeries)
with open(outfilename, 'w') as outfile:
    json.dump(data.data, outfile, indent=2)
