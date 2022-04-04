import json

class jsonTestData:
    def __init__(self, database):
        self.database = database

        self.mstdbReferences = []

        # Set data filters empty
        self.seriesStatusFilter = []
        self.seriesTypeFilter = []
        self.seriesElementsFilter = []
        self.sampleStatusFilter = []

        # Call data filter tool (with empty filters)
        self.filter()

    def filter(self):
        f = open(self.database,)
        try:
            self.data = json.load(f)
            f.close()
        except:
            print('Failed to load test file')
            exit()
        self.nRefs = 0
        self.nSources = 0
        self.nSeries = 0
        self.nSamples = 0

        delSeriesList = []
        delSourceList = []

        for sourceName in self.data['sources'].keys():
            source = self.data['sources'][sourceName]
            sourceSeries  = 0
            sourceSamples = 0
            mstdbParentRefs = source['mstdb references']
            for testType in source['tests'].keys():
                if testType in ['vapor pressures','solubility limits','heat capacities']:
                    for seriesName in source['tests'][testType].keys():
                        series = source['tests'][testType][seriesName]
                        seriesSamples = 0
                        # Calculate status of series before deciding if samples are included
                        seriesStatus = 'unknown'
                        for sample in series['samples'].keys():
                            # If no status, write incomplete status
                            try:
                                series['samples'][sample]['status']
                            except KeyError:
                                series['samples'][sample]['status'] = 'incomplete'
                            sampleStatus = series['samples'][sample]['status']
                            if sampleStatus == 'pass':
                                if seriesStatus == 'unknown':
                                    seriesStatus = 'pass'
                                elif seriesStatus == 'fail':
                                    seriesStatus = 'partial'
                            elif sampleStatus == 'fail':
                                if seriesStatus == 'unknown':
                                    seriesStatus = 'fail'
                                elif seriesStatus == 'pass':
                                    seriesStatus = 'partial'
                            else:
                                seriesStatus = 'incomplete'
                                sampleStatus = 'incomplete'
                            if (sampleStatus in self.sampleStatusFilter) or (not self.sampleStatusFilter):
                                seriesSamples += 1
                        series['series status'] = seriesStatus
                        series['series type'] = testType

                        # Check series filters for series inclusion
                        includeSeries = True
                        if (seriesStatus not in self.seriesStatusFilter) and self.seriesStatusFilter:
                            includeSeries = False
                        if (testType not in self.seriesTypeFilter) and self.seriesTypeFilter:
                            includeSeries = False
                        if seriesSamples <= 0:
                            includeSeries = False
                        for element in self.seriesElementsFilter:
                            if element not in series['composition'].keys():
                                includeSeries = False
                        # If series included, update source variables
                        if includeSeries:
                            sourceSeries  += 1
                            sourceSamples += seriesSamples
                        else:
                            delSeriesList.append([sourceName,testType,seriesName])
                elif testType == 'phase transitions':
                    for seriesName in source['tests'][testType].keys():
                        series = source['tests'][testType][seriesName]
                        # For phase transitions, series and samples are equivalent
                        seriesSamples = 0
                        seriesStatus = 'incomplete'
                        if 'status' not in list(series.keys()):
                            seriesStatus = 'incomplete'
                        elif series['status'] == 'pass':
                            seriesStatus = 'pass'
                        elif series['status'] == 'fail':
                            seriesStatus = 'fail'
                        # Apply the sample filters to series for phase transitions
                        if (series['status'] in self.sampleStatusFilter) or (not self.sampleStatusFilter):
                            seriesSamples += 1
                        series['series status'] = seriesStatus

                        # Check series filters for series inclusion
                        includeSeries = True
                        if (seriesStatus not in self.seriesStatusFilter) and self.seriesStatusFilter:
                            includeSeries = False
                        if (testType not in self.seriesTypeFilter) and self.seriesTypeFilter:
                            includeSeries = False
                        if seriesSamples <= 0:
                            includeSeries = False
                        for element in self.seriesElementsFilter:
                            if element not in series['composition'].keys():
                                includeSeries = False
                        # If series included, update source variables
                        if includeSeries:
                            sourceSeries  += 1
                            sourceSamples += seriesSamples
                        else:
                            delSeriesList.append([sourceName,testType,seriesName])

            # Check source filters for source inclusion
            includeSource = True
            if sourceSeries <= 0:
                includeSource = False
            if sourceSamples <= 0:
                includeSource = False
            # If source included, update total variables
            if includeSource:
                self.nSources += 1
                self.nSeries  += sourceSeries
                self.nSamples += sourceSamples
            else:
                delSourceList.append(sourceName)

        for item in delSeriesList:
            del self.data['sources'][item[0]]['tests'][item[1]][item[2]]
        for item in delSourceList:
            del self.data['sources'][item]

        self.nRefs = len(self.mstdbReferences)
        print(self.nSources)
        print(self.nSeries)
        print(self.nSamples)
