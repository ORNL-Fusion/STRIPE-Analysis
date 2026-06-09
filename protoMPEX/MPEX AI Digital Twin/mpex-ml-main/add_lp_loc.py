#-------------------------------------------------------------------------------
#  add_lp_loc.py
#  Python script to add langmuir probe location data to the SQL database.
#-------------------------------------------------------------------------------

import MDSplus
import psycopg2
import argparse
import math

#-------------------------------------------------------------------------------
#  Run main program if the script is run directly.
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    command_line_parser = argparse.ArgumentParser()

    command_line_parser.add_argument('-m',
                                     '--mdsplusserver',
                                     action='store',
                                     default='mpexserver.ornl.gov',
                                     dest='mdsplusserver',
                                     help='IP address or URL of the mdsplus server',
                                     metavar='MODULE_NAME')

    args = vars(command_line_parser.parse_args())

#  Remove empty arguments
    for key in [key for key in args if args[key] == None]:
        del args[key]

    connection = MDSplus.Connection(args['mdsplusserver'])

    sql = psycopg2.connect(dbname='2ya', port=5432, host='localhost', user='2ya')

    with sql:
        with sql.cursor() as cursor:
            cursor.execute('SELECT shotnumber FROM lp;')

            data = []

            for record in cursor:
                print(record[0])
                try:
                    connection.get('TreeOpen("mpex",$)', record[0])

                    data.append({'shotnumber' : record[0],
                                 'axial_loc'  : float(connection.get('ANALYZED.DLP.SETUP:AXIAL_LOC')),
                                 'rad_loc'    : float(connection.get('ANALYZED.DLP.SETUP:RAD_LOC'))})
                except Exception as e:
                    None

            for record in data:
                print(record, not math.isnan(record['axial_loc']) and not math.isnan(record['rad_loc']))
                if not math.isnan(record['axial_loc']) and not math.isnan(record['rad_loc']):
                    cursor.execute('UPDATE lp SET axial_loc = {axial_loc}, rad_loc = {rad_loc} WHERE shotnumber = {shotnumber}'.format(**record))
