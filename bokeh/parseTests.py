import json

class jsonTestData:
    def __init__(self, database):
        f = open(database,)
        try:
            self.data = json.load(f)
            f.close()
        except:
            print('Failed to load test file')
            exit()

        self.mstdbReferences = []
        self.experimentalReferences = []
        self.testSeries = []

        for source in self.data['sources'].keys():
            mstdbParentRefs = self.data['sources'][source]['mstdb references']
            for testType in self.data['sources'][source]['tests'].keys():
                if testType in ['vapor pressures','solubility limits','heat capacities']:
                    for series in self.data['sources'][source]['tests'][testType].keys():
                        # Calculate status of series before deciding if samples are included
                        seriesStatus = 'unknown'
                        for sample in self.data['sources'][source]['tests'][testType][series]['samples'].keys():
                            if 'status' not in list(self.data['sources'][source]['tests'][testType][series]['samples'][sample].keys()):
                                seriesStatus = 'incomplete'
                                break
                            if self.data['sources'][source]['tests'][testType][series]['samples'][sample]['status'] == 'pass':
                                if seriesStatus == 'unknown':
                                    seriesStatus = 'pass'
                                elif seriesStatus == 'fail':
                                    seriesStatus = 'partial'
                            elif self.data['sources'][source]['tests'][testType][series]['samples'][sample]['status'] == 'fail':
                                if seriesStatus == 'unknown':
                                    seriesStatus = 'fail'
                                elif seriesStatus == 'pass':
                                    seriesStatus = 'partial'
                            else:
                                seriesStatus = 'incomplete'
                                break
                        self.data['sources'][source]['tests'][testType][series]['series status'] = seriesStatus
                        self.data['sources'][source]['tests'][testType][series]['series type'] = testType
                        # For included tests, references to respective lists
                        for ref in mstdbParentRefs:
                            if ref not in self.mstdbReferences:
                                self.mstdbReferences.append(ref)
                        if self.data['sources'][source] not in self.experimentalReferences:
                            self.experimentalReferences.append(self.data['sources'][source])

                        # Add test to series list
                        self.testSeries.append(dict([(series,self.data['sources'][source]['tests'][testType][series])]))
                elif testType == 'phase transitions':
                    for series in self.data['sources'][source]['tests'][testType].keys():
                        seriesStatus = 'incomplete'
                        if 'status' not in list(self.data['sources'][source]['tests'][testType][series].keys()):
                            seriesStatus = 'incomplete'
                        elif self.data['sources'][source]['tests'][testType][series]['status'] == 'pass':
                            seriesStatus = 'pass'
                        elif self.data['sources'][source]['tests'][testType][series]['status'] == 'fail':
                            seriesStatus = 'fail'
                        self.data['sources'][source]['tests'][testType][series]['series status'] = seriesStatus
                        # For included tests, references to respective lists
                        for ref in mstdbParentRefs:
                            if ref not in self.mstdbReferences:
                                self.mstdbReferences.append(ref)
                        if self.data['sources'][source] not in self.experimentalReferences:
                            self.experimentalReferences.append(self.data['sources'][source])

                        # Add test to series list
                        self.testSeries.append(dict([(series,self.data['sources'][source]['tests'][testType][series])]))
