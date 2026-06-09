#-------------------------------------------------------------------------------
#  data_avail.py
#  Python program to populate and query the database of data avialability.
#-------------------------------------------------------------------------------

import argparse
import json
import psycopg2
import shot_data
import MDSplus

#-------------------------------------------------------------------------------
#  Add rows to database.
#-------------------------------------------------------------------------------
def populate(connection, config, start_shot, end_shot):
    preamble = 'REPLACE INTO avail (shotnumber'

    sql = psycopg2.connect(dbname='2ya', port=5432, host='localhost')

    for i in range(args['start_shot'], args['end_shot']):
        shot = shot_data.mpex_shot(connection, i, config)
        shot.close()

        channel_names = shot.get_channel_names()

        sqlcommand = preamble
        middle = ') VALUES ({}'.format(i)

        for name in channel_names:
            preamble += ',{}'.format(name)
            middle += ',TRUE'

        with sql:
            with sql.cursor() as cursor:
                cursor.execute('{}{});'.format(preamble, middle))

    sql.close()

#-------------------------------------------------------------------------------
#  Add rows to database.
#-------------------------------------------------------------------------------
def query(connection):
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
    command_line_parser.add_argument('-q'
                                     '--query',
                                     action='store',
                                     dest='query',
                                     help='SQL Query command.',
                                     metavar='QUERY')
    command_line_parser.add_argument('-s',
                                     '--start_shot',
                                     action='store',
                                     default=1,
                                     type=int,
                                     dest='start_shot',
                                     help='Starting shot number to populate database.',
                                     metavar='START_SHOT')
    command_line_parser.add_argument('-e',
                                     '--end_shot',
                                     action='store',
                                     default=1,
                                     type=int,
                                     dest='end_shot',
                                     help='Ending shot number to populate database.',
                                     metavar='END_SHOT')

    args = vars(command_line_parser.parse_args())

#  Remove empty arguments
    for key in [key for key in args if args[key] == None]:
        del args[key]

    with open(args['config'], 'r') as json_ref:
        config = json.load(json_ref)

    connection = MDSplus.Connection(args['mdsplusserver'])

    if 'start_shot' in args and 'end_shot' in args:
        populate(connection, config, args['start_shot'], args['end_shot'])
    elif 'query' in args:
        query(connection)
    else:
        print('Unknown task:')
        print('  set start_shot and end_shot or query in commandline arguments.')
