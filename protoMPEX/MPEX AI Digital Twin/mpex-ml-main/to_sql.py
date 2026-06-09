#-------------------------------------------------------------------------------
#  to_sql.py
#  Python class to convert shot data to SQL
#-------------------------------------------------------------------------------

import shot_data
import MDSplus
import psycopg2
import argparse
import json
import numpy

#-------------------------------------------------------------------------------
#  Format an numpy array into an SLQ array.
#
#  param[in] array Array to convert.
#-------------------------------------------------------------------------------
def to_sql_array(array):
    command = '{}'.format(array[0])
    for value in array[1:]:
        command += ', {}'.format(value)
    return command

#-------------------------------------------------------------------------------
#  Format an numpy array into an SLQ array.
#
#  param[in] array Array to convert.
#-------------------------------------------------------------------------------
def to_sql_2darray(array):
    command = 'ARRAY[{}'.format(to_sql_array(array[0]))
    for value in array[1:]:
        command += ', {}'.format(to_sql_array(value))
    command += ']'
    return command

#-------------------------------------------------------------------------------
#  Class to interface with SQL Database
#-------------------------------------------------------------------------------
class mds_to_sql:
#-------------------------------------------------------------------------------
#  Initialize a mds_to_sql instance.
#
#  This loads channel data and time base for a channel.
#
#  param[inout] self       A mds_to_sql instance.
#  param[in]    connection A connection to the MDSPlus server.
#  param[in]    number     Shot number for the data.
#  param[in]    config     Channel configureation.
#-------------------------------------------------------------------------------
    def __init__(self, connection, shot_number, config):
        self.shot = shot_data.mpex_shot(connection, shot_number, config)
        self.shot.close()

        self.sql = psycopg2.connect(dbname='2ya', port=5432, host='localhost')

#-------------------------------------------------------------------------------
#  Finalize a mds_to_sql instance.
#
#  param[inout] self A mds_to_sql instance.
#-------------------------------------------------------------------------------
    def __del__(self):
        self.sql.close()

#-------------------------------------------------------------------------------
#  Insert timeslice.
#
#  param[inout] self A mds_to_sql instance.
#-------------------------------------------------------------------------------
    def insert(self, start_time, end_time):
        try:
            preamble = 'INSERT INTO mpex (shotnumber, timeslice'
            middle = ') VALUES ('
            ending = ');'

            timeslice = int(1000.0*(end_time + start_time)/2.0)

            slice = shot_data.mpex_time_slice(self.shot, start_time, end_time)

            command = preamble
            for group in slice.channels:
                for channel in slice.channels[group]:
                    command += ', {}'.format(channel)
            command += '{}{}, {}'.format(middle, slice.shot_number, timeslice)
            for group in slice.channels:
                for channel in slice.channels[group]:
                    if channel == 'B_TOTAL':
                        command += ', {}'.format(to_sql_2darray(slice.channels[group][channel]))
                    else:
                        command += ', {}'.format(slice.channels[group][channel])
            command += ending

            with self.sql:
                with self.sql.cursor() as cursor:
                    cursor.execute(command)

            print(slice.shot_number, timeslice, start_time, end_time)
        except Exception as e:
            None

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
    command_line_parser.add_argument('-c',
                                     '--config',
                                     action='store',
                                     required=True,
                                     dest='config',
                                     help='mdsplus channel configuation',
                                     metavar='CONFIG')
    command_line_parser.add_argument('-s',
                                     '--startshotnumber',
                                     action='store',
                                     default=22167,
                                     type=int,
                                     dest='startshotnumber',
                                     help='Start shot number to grab mdsplus data from.',
                                     metavar='STARTSHOTNUMBER')
    command_line_parser.add_argument('-e',
                                     '--endshotnumber',
                                     action='store',
                                     default=25414,
                                     type=int,
                                     dest='endshotnumber',
                                     help='End shot number to grab mdsplus data from.',
                                     metavar='ENDSHOTNUMBER')
    command_line_parser.add_argument('-st',
                                     '--starttime',
                                     action='store',
                                     default=4.1,
                                     type=float,
                                     dest='starttime',
                                     help='Starting time to slice data.',
                                     metavar='STARTTIME')
    command_line_parser.add_argument('-et',
                                     '--endtime',
                                     action='store',
                                     default=5.0,
                                     type=float,
                                     dest='endtime',
                                     help='Ending time to slice data.',
                                     metavar='ENDTIME')
    command_line_parser.add_argument('-nt',
                                     '--numtimes',
                                     action='store',
                                     default=10,
                                     type=int,
                                     dest='numtimes',
                                     help='Number of slices to make.',
                                     metavar='NUMTIMES')

    args = vars(command_line_parser.parse_args())

#  Remove empty arguments
    for key in [key for key in args if args[key] == None]:
        del args[key]

    with open(args['config'], 'r') as json_ref:
        config = json.load(json_ref)

    connection = MDSplus.Connection(args['mdsplusserver'])

    for shotnumber in range(args['startshotnumber'], args['endshotnumber'] + 1):
        sql = mds_to_sql(connection, shotnumber, config)

        times = numpy.linspace(args['starttime'], args['endtime'], args['numtimes'])
        for i in range(args['numtimes'] - 1):
            sql.insert(times[i], times[i + 1]);
