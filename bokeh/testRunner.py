import json
import parseTests
import subprocess
import copy
import collections.abc

def runVaporPressures(series,name):
    inputScript = f'vaporPressures/{name}.ti'
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
    # Check for alternative units
    if "temperature unit" in series.keys():
        tseriesunit = series["temperature unit"]
    else:
        tseriesunit = tunit
    vaporPressureScale = 1
    if "vapor pressure unit" in series.keys():
        if series["vapor pressure unit"] == 'mmHG':
            vaporPressureScale = 1/760
    # Write input script
    with open(inputScript, 'w') as inputFile:
        inputFile.write('! Python-generated input file for Thermochimica\n')
        inputFile.write(f'data file         = {datapath}\n')
        inputFile.write(f'temperature unit  = {tseriesunit}\n')
        inputFile.write(f'pressure unit     = {punit}\n')
        inputFile.write(f'mass unit         = {munit}\n')
        inputFile.write(f'nEl               = {nElements}\n')
        inputFile.write(f'iEl               = {" ".join([str(atomic_number_map.index(element)+1) for element in elements])}\n')
        inputFile.write(f'nCalc             = {len(list(series["samples"].keys()))}\n')
        for sample in series['samples']:
            s = series['samples'][sample]
            inputFile.write(f'{s["temperature"]} {press} {compString}\n')
    # Run calculations
    subprocess.run([thermopath,inputScript])
    # Get output from Thermochimica
    try:
        f = open(outjsonpath,)
        out = json.load(f)
        f.close()
    except:
        print('Failed to open Thermochimica output file')
        for sample in series['samples']:
            series['samples'][sample]['status'] = 'fail'
        return
    # Check output against experiment
    for sample in series['samples']:
        s = series['samples'][sample]
        try:
            o = out[sample]['solution phases']['gas_ideal']['species']
        except KeyError:
            s['status'] = 'fail'
            continue
        # Check list of excluded species
        adjustedTotal = 1
        try:
            for species in series['exclude']:
                adjustedTotal -= o[species]['mole fraction']*press
        except:
            pass
        sampleStatus = 'pass'
        for species in s['partial pressures']:
            try:
                # Calculate bounds
                lb = (10**(-logErr)) * s['partial pressures'][species] * vaporPressureScale
                ub = (10**( logErr)) * s['partial pressures'][species] * vaporPressureScale
                calculated = o[species]['mole fraction']*press / adjustedTotal
                if not (calculated >= lb and calculated <= ub):
                    sampleStatus = 'fail'
            except KeyError:
                sampleStatus = 'fail'
        s['status'] = sampleStatus

def runPhaseTransitions(series,name):
    inputScript = f'phaseTransitions/{name}.ti'
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
    target = series['target temperature']
    err = series['error']
    plo = series['phases low']
    phi = series['phases high']
    phaseMolTol = 1e-6
    # Check for alternative units
    if "temperature unit" in series.keys():
        tseriesunit = series["temperature unit"]
    else:
        tseriesunit = tunit
    # Write input script
    with open(inputScript, 'w') as inputFile:
        inputFile.write('! Python-generated input file for Thermochimica\n')
        inputFile.write(f'data file         = {datapath}\n')
        inputFile.write(f'temperature unit  = {tseriesunit}\n')
        inputFile.write(f'pressure unit     = {punit}\n')
        inputFile.write(f'mass unit         = {munit}\n')
        inputFile.write(f'nEl               = {nElements}\n')
        inputFile.write(f'iEl               = {" ".join([str(atomic_number_map.index(element)+1) for element in elements])}\n')
        inputFile.write(f'nCalc             = 2\n')
        # Run once at low temperature limit and once at high
        inputFile.write(f'{target - err} {press} {compString}\n')
        inputFile.write(f'{target + err} {press} {compString}\n')
    # Run calculations
    subprocess.run([thermopath,inputScript])
    # Get output from Thermochimica
    try:
        f = open(outjsonpath,)
        out = json.load(f)
        f.close()
    except:
        print('Failed to open Thermochimica output file')
        series['status'] = 'fail'
        return
    # Check output against experiment
    # Low temperature check
    series['status'] = 'pass'
    # Check if run converged
    try:
        out["1"]["solution phases"].keys()
    except KeyError:
        print('Low temperature run did not converge')
        series['status'] = 'fail'
        return
    for phase in out["1"]["solution phases"].keys():
        p = out["1"]["solution phases"][phase]
        # If phase is supposed to be there, make sure it has non-zero moles
        if phase in plo:
            if p['moles'] <= 0.0:
                series['status'] = 'fail'
        # Otherwise, should have zero
        else:
            if p['moles'] > phaseMolTol:
                series['status'] = 'fail'
    for phase in out["1"]["pure condensed phases"].keys():
        p = out["1"]["pure condensed phases"][phase]
        # If phase is supposed to be there, make sure it has non-zero moles
        if phase in plo:
            if p['moles'] <= 0.0:
                series['status'] = 'fail'
        # Otherwise, should have zero
        else:
            if p['moles'] > phaseMolTol:
                series['status'] = 'fail'
    # Now reverse: if it is supposed to be there, ensure it is
    for phase in plo:
        try:
            if out["1"]["solution phases"][phase]['moles'] <= 0.0:
                series['status'] = 'fail'
        except KeyError:
            # If not in solution phases, check pure condensed
            try:
                if out["1"]["pure condensed phases"][phase]['moles'] <= 0.0:
                    series['status'] = 'fail'
            except KeyError:
                print(f'Phase {phase} not found in Thermochimica output')
                series['status'] = 'fail'

    # High temperature check
    # Check if run converged
    try:
        out["2"]["solution phases"].keys()
    except KeyError:
        print('High temperature run did not converge')
        series['status'] = 'fail'
        return
    for phase in out["2"]["solution phases"].keys():
        p = out["2"]["solution phases"][phase]
        # If phase is supposed to be there, make sure it has non-zero moles
        if phase in phi:
            if p['moles'] <= 0.0:
                series['status'] = 'fail'
        # Otherwise, should have zero
        else:
            if p['moles'] > phaseMolTol:
                series['status'] = 'fail'
    for phase in out["2"]["pure condensed phases"].keys():
        p = out["2"]["pure condensed phases"][phase]
        # If phase is supposed to be there, make sure it has non-zero moles
        if phase in phi:
            if p['moles'] <= 0.0:
                series['status'] = 'fail'
        # Otherwise, should have zero
        else:
            if p['moles'] > phaseMolTol:
                series['status'] = 'fail'
    # Now reverse: if it is supposed to be there, ensure it is
    for phase in phi:
        try:
            if out["2"]["solution phases"][phase]['moles'] <= 0.0:
                series['status'] = 'fail'
        except KeyError:
            # If not in solution phases, check pure condensed
            try:
                if out["2"]["pure condensed phases"][phase]['moles'] <= 0.0:
                    series['status'] = 'fail'
            except KeyError:
                print(f'Phase {phase} not found in Thermochimica output')
                series['status'] = 'fail'

def runSolubilityLimits(series,name):
    inputScript = f'solubilityLimits/{name}.ti'
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
    err = series['error']
    phase = series['phase']
    fracType = series['fraction type']
    # Check for alternative units
    if "temperature unit" in series.keys():
        tseriesunit = series["temperature unit"]
    else:
        tseriesunit = tunit
    # Write input script
    with open(inputScript, 'w') as inputFile:
        inputFile.write('! Python-generated input file for Thermochimica\n')
        inputFile.write(f'data file         = {datapath}\n')
        inputFile.write(f'temperature unit  = {tseriesunit}\n')
        inputFile.write(f'pressure unit     = {punit}\n')
        inputFile.write(f'mass unit         = {munit}\n')
        inputFile.write(f'nEl               = {nElements}\n')
        inputFile.write(f'iEl               = {" ".join([str(atomic_number_map.index(element)+1) for element in elements])}\n')
        inputFile.write(f'nCalc             = {len(list(series["samples"].keys()))}\n')
        for sample in series['samples']:
            s = series['samples'][sample]
            inputFile.write(f'{s["temperature"]} {press} {compString}\n')
    # Run calculations
    subprocess.run([thermopath,inputScript])
    # Get output from Thermochimica
    try:
        f = open(outjsonpath,)
        out = json.load(f)
        f.close()
    except:
        print('Failed to open Thermochimica output file')
        for sample in series['samples']:
            series['samples'][sample]['status'] = 'fail'
        return
    # Check output against experiment
    for sample in series['samples']:
        s = series['samples'][sample]
        try:
            o = out[sample]['solution phases'][phase][fracType]
        except KeyError:
            s['status'] = 'fail'
            continue
        s['status'] = 'pass'
        # Make sure phase is stable
        if out[sample]['solution phases'][phase]['moles'] <= 0.0:
            s['status'] = 'fail'
        # Calculate adjusted total pairs for multiple coordination cases
        multipleCoordPairs = ['Al2Cl6','Pu2Cl6','Be2F4']
        totalPairs = 1
        for pair in multipleCoordPairs:
            if pair in o.keys():
                totalPairs += o[pair]['mole fraction']
        for species in s['fractions']:
            try:
                # Calculate bounds
                lb = s['fractions'][species] - err
                ub = s['fractions'][species] + err
                calculated = o[species]['mole fraction']
                # Need to handle multiple coordination number cases specially
                if species == 'AlCl3':
                    calculated = calculated + 2*o['Al2Cl6']['mole fraction']
                if species == 'PuCl3':
                    calculated = calculated + 2*o['Pu2Cl6']['mole fraction']
                if species == 'BeF2':
                    calculated = calculated + 2*o['Be2F4']['mole fraction']
                calculated = calculated / totalPairs
                if not (calculated >= lb and calculated <= ub):
                    s['status'] = 'fail'
            except KeyError:
                s['status'] = 'fail'

def runHeatCapacities(series,name):
    inputScript = f'heatCapacities/{name}.ti'
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
    relErr = series['relative error']
    # Check for alternative units
    if "temperature unit" in series.keys():
        tseriesunit = series["temperature unit"]
    else:
        tseriesunit = tunit
    # Write input script
    with open(inputScript, 'w') as inputFile:
        inputFile.write('! Python-generated input file for Thermochimica\n')
        inputFile.write(f'heat capacity     = .TRUE.\n')
        inputFile.write(f'data file         = {datapath}\n')
        inputFile.write(f'temperature unit  = {tseriesunit}\n')
        inputFile.write(f'pressure unit     = {punit}\n')
        inputFile.write(f'mass unit         = {munit}\n')
        inputFile.write(f'nEl               = {nElements}\n')
        inputFile.write(f'iEl               = {" ".join([str(atomic_number_map.index(element)+1) for element in elements])}\n')
        inputFile.write(f'nCalc             = {len(list(series["samples"].keys()))}\n')
        for sample in series['samples']:
            s = series['samples'][sample]
            inputFile.write(f'{s["temperature"]} {press} {compString}\n')
    # Run calculations
    subprocess.run([thermopath,inputScript])
    # Get output from Thermochimica
    try:
        f = open(outjsonpath,)
        out = json.load(f)
        f.close()
    except:
        print('Failed to open Thermochimica output file')
        for sample in series['samples']:
            series['samples'][sample]['status'] = 'fail'
        return
    # Check output against experiment
    for sample in series['samples']:
        s = series['samples'][sample]
        try:
            o = out[sample]
        except KeyError:
            s['status'] = 'fail'
            continue
        sampleStatus = 'pass'
        try:
            # Calculate bounds
            lb = (1 - relErr) * s['heat capacity']
            ub = (1 + relErr) * s['heat capacity']
            calculated = o['heat capacity']
            if not (calculated >= lb and calculated <= ub):
                sampleStatus = 'fail'
        except KeyError:
            sampleStatus = 'fail'
        s['status'] = sampleStatus

# Set file names for input/output
infilename  = 'validationData.json'
outfilename = 'validationData.json'

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

# utility to recursively update a dict with another dict
def updateDictWithDict(d, u):
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = updateDictWithDict(d.get(k, {}), v)
        else:
            d[k] = v
    return d

def run(data):
    for sourceName in data.data['sources']:
        source = data.data['sources'][sourceName]
        # Loop over all associated test series and execture
        for testType in source['tests']:
            for series in source['tests'][testType]:
                print(series)
                name = f'{sourceName}-{series}'
                currentSeries = source['tests'][testType][series]
                if 'enabled' in currentSeries.keys():
                    if currentSeries['enabled'] in ['false', 'False', 'FALSE', 0]:
                        print('disabled')
                        continue
                if testType in ['vapor pressures']:
                    runVaporPressures(currentSeries,name)
                    for sample in currentSeries['samples']:
                        print(f"{sample}: {currentSeries['samples'][sample]['status']}")
                elif testType == 'phase transitions':
                    runPhaseTransitions(currentSeries,name)
                    print(currentSeries['status'])
                elif testType in ['solubility limits']:
                    runSolubilityLimits(currentSeries,name)
                    for sample in currentSeries['samples']:
                        print(f"{sample}: {currentSeries['samples'][sample]['status']}")
                elif testType in ['heat capacities']:
                    runHeatCapacities(currentSeries,name)
                    for sample in currentSeries['samples']:
                        print(f"{sample}: {currentSeries['samples'][sample]['status']}")

def runAll():
    # Get data
    data = parseTests.jsonTestData(infilename)

    run(data)

    with open(outfilename, 'w') as outfile:
        json.dump(data.data, outfile, indent=2)

def runNew():
    # Get data
    data = parseTests.jsonTestData(infilename)
    fullData = copy.deepcopy(data)
    data.seriesStatusFilter = ['incomplete']
    data.filter()

    run(data)

    updateDictWithDict(fullData.data, data.data)

    with open(outfilename, 'w') as outfile:
        json.dump(fullData.data, outfile, indent=2)

def runSource(sourceNames):
    # sourceNames is array of source names to include
    # e.g. ['vorobei-1973','thoma-1959']
    # Get data
    data = parseTests.jsonTestData(infilename)
    fullData = copy.deepcopy(data)
    data.sourceFilter = sourceNames
    data.filter()

    run(data)

    updateDictWithDict(fullData.data, data.data)

    with open(outfilename, 'w') as outfile:
        json.dump(fullData.data, outfile, indent=2)

def runTestType(testTypes):
    # testTypes is array of type names to include
    # e.g. ['phase transitions','heat capacities']
    # Set in/out-put names
    infilename  = 'validationData.json'
    outfilename = 'phaseTransitions.json'
    # Get data
    data = parseTests.jsonTestData(infilename)
    fullData = copy.deepcopy(data)
    data.seriesTypeFilter = testTypes
    data.filter()

    run(data)

    # updateDictWithDict(fullData.data, data.data)

    with open(outfilename, 'w') as outfile:
        json.dump(fullData.data, outfile, indent=2)


runTestType(['phase transitions'])
