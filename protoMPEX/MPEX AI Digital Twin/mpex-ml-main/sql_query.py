import psycopg2
import tensorflow
import json
import numpy as np
import scipy.io as sio

def transpose_data(data):
    tensor = []
    for row in data:
        temprow = []
        for key in row:
            if type(row[key]) is list:
                temprow += row[key]
            else:
                temprow.append(row[key])
        tensor.append(tensorflow.constant([temprow]))

    return tensorflow.concat(tensor, 0)

def to_json(name, data):
    transpose_data = {}
    for row in data:
        for key in row:
            if key in transpose_data:
                transpose_data[key].append(row[key])
            else:
                transpose_data[key] = [row[key]]
    with open('{}.json'.format(name), 'w') as json_ref:
        json.dump(transpose_data, json_ref, indent=4)

sql = psycopg2.connect(dbname='2ya', port=5432, host='localhost', user='2ya')

indata = []
outdata = []

with sql:
    with sql.cursor() as cursor:
        cursor.execute('SELECT mpex.pwr_hel, mpex.pwr_28, mpex.pwr_icrf, mpex.b_total, mpex.pg2_pres, mpex.pg3_pres, mpex.pg4_pres, lp.axial_loc, lp.rad_loc FROM mpex JOIN lp ON mpex.shotnumber = lp.shotnumber JOIN bad ON mpex.shotnumber != bad.shotnumber WHERE mpex.heat_peak_ce IS NOT NULL AND mpex.heat_peak IS NOT NULL AND mpex.b_total IS NOT NULL AND mpex.pg2_pres IS NOT NULL AND mpex.pg3_pres IS NOT NULL AND mpex.pg4_pres IS NOT NULL AND mpex.kte IS NOT NULL AND mpex.ne IS NOT NULL AND lp.axial_loc IS NOT NULL AND lp.rad_loc IS NOT NULL;')

        for record in cursor:
            indata.append({'pwr_hel'   : record[0],
                           'pwr_28'    : record[1],
                           'pwr_icrf'  : record[2],
                           'b_total'   : record[3],
                           'pg2_pres'  : record[4],
                           'pg3_pres'  : record[5],
                           'pg4_pres'  : record[6],
                           'axial_loc' : record[7],
                           'rad_loc'   : record[8]
                          })

        cursor.execute('SELECT mpex.heat_peak_ce, mpex.heat_peak, mpex.kte, mpex.ne FROM mpex JOIN lp ON mpex.shotnumber = lp.shotnumber JOIN bad ON mpex.shotnumber != bad.shotnumber WHERE mpex.heat_peak_ce IS NOT NULL AND mpex.heat_peak IS NOT NULL AND mpex.b_total IS NOT NULL AND mpex.pg2_pres IS NOT NULL AND mpex.pg3_pres IS NOT NULL AND mpex.pg4_pres IS NOT NULL AND mpex.kte IS NOT NULL AND mpex.ne IS NOT NULL AND lp.axial_loc IS NOT NULL AND lp.rad_loc IS NOT NULL;')

        for record in cursor:
            outdata.append({'heat_peak_ce' : record[0],
                            'heat_peak'    : record[1],
                            'kte'          : record[2],
                            'ne'           : record[3]
                           })

to_json('indata', indata)
to_json('outdata', outdata)

outdata = transpose_data(outdata)
indata = transpose_data(indata)
print(tensorflow.shape(indata))
print(tensorflow.shape(outdata))


indat = np.array(indata)
outdat = np.array(outdata)


sio.savemat('indat.mat', {'indat': indat})
sio.savemat('outdat.mat', {'outdat': outdat})
