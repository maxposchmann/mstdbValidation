import json
import csv

counter = 0
jsondata = dict([])
with open('../references/romberger-1972-f5-ii.csv', newline='') as csvfile:
    csvdata = csv.reader(csvfile, delimiter=',')
    for row in csvdata:
        counter += 1
        # heat capacities
        # jsondata[str(counter)] = dict([("temperature",float(row[0])),("heat capacity",float(row[1]))])
        # solubility limits
        species = 'BeF2'
        jsondata[str(counter)] = dict([("temperature",float(row[1])),("fractions",dict([(species,float(row[0]))]))])

with open('out.json', 'w') as outfile:
    json.dump(jsondata, outfile, indent=2)
