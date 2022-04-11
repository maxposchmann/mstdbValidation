import json
import csv

counter = 0
jsondata = dict([])
with open('../references/gaune-escard-1999-f1.csv', newline='') as csvfile:
    csvdata = csv.reader(csvfile, delimiter=',')
    for row in csvdata:
        counter += 1
        jsondata[str(counter)] = dict([("temperature",float(row[0])),("heat capacity",float(row[1]))])

with open('out.json', 'w') as outfile:
    json.dump(jsondata, outfile, indent=2)
