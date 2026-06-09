#-------------------------------------------------------------------------------
#  shot_data.py
#  Python class to represent the data for single shot.
#-------------------------------------------------------------------------------

import MDSplus
import matplotlib.pyplot
import numpy
import argparse
import json
import copy

#-------------------------------------------------------------------------------
#  Class to represent a single channel.
#-------------------------------------------------------------------------------
class mpex_channel:
#-------------------------------------------------------------------------------
#  Initialize a mpex_channel.
#
#  This loads channel data and time base for a channel.
#
#  param[inout] self       A mpex_channel instance.
#  param[in]    connection A connection to the MDSPlus server.
#  param[in]    name       Tree name of the channel.
#  param[in]    operations Operations to transform the data.
#-------------------------------------------------------------------------------
    def __init__(self, connection, name, operations):
        self.value = connection.get(name).data()
        self.name = name

        if operations != None:
            for operation in operations:
                if operation == 'flatten':
                    self.value = numpy.ndarray.flatten(self.value)
                if operation == 'transpose':
                    self.value = numpy.transpose(self.value)

#-------------------------------------------------------------------------------
#  Plot a mpex_channel:
#
#  Mpex channels can be 1, 2, or 3 dimensional + time.
#
#  param[in] self     A mpex_channel instance.
#  param[in] axis     Area to plot data to.
#  param[in] timebase Time axis.
#  param[in] position Position axis.
#  param[in] min_time Minimum plot time.
#  param[in] max_time Maximum plot time.
#-------------------------------------------------------------------------------
    def plot(self, axis, timebase, position, min_time, max_time):
        axis.set_title(self.name)
        axis.set_xlim([min_time, max_time])

        if timebase != None:
            if position != None:
                axis.imshow(numpy.transpose(self.value[:,1,:]), extent=(numpy.ndarray.min(timebase.value),numpy.ndarray.max(timebase.value),numpy.ndarray.min(position.value),numpy.ndarray.max(position.value)), aspect='auto')
            else:
                try:
                    axis.plot(timebase.value, self.value)
                except ValueError:
                    axis.plot(timebase.value, timebase.value)
                    print('{} Failed due to dimensionality mismatch.'.format(self.name))
        else:
            print('{} Failed due to missing timebase.'.format(self.name))

#-------------------------------------------------------------------------------
#  Print a mpex_channel.
#
#  param[in] self A mpex_channel instance.
#-------------------------------------------------------------------------------
    def print(self):
        print(self.name)

#-------------------------------------------------------------------------------
#  Slice a mpex_channel.
#
#  Slices work by averaging data with in the time interval. The time window is
#  too small for the time base, the average is taken from the first ajacent
#  points before and after the time interval.
#
#  param[in] self       A mpex_channel instance.
#  param[in] start_time Start time of the slice.
#  param[in] end_time   End time of the slice.
#  param[in] timebase   Timebase of the channel.
#-------------------------------------------------------------------------------
    def slice(self, start_time, end_time, timebase):
        time_mask = numpy.logical_and(timebase >= start_time, timebase <= end_time)
        temp_value = self.value[time_mask]

        if temp_value.size > 0:
            temp = numpy.sum(temp_value, 0)/numpy.size(timebase[time_mask])
            return numpy.sum(temp_value, 0)/numpy.size(timebase[time_mask])
        else:
            lower = numpy.where(timebase < start_time)[-1][-1]
            upper = numpy.where(timebase > end_time)[0][0]

            slope = (self.value[upper] - self.value[lower])/(timebase[upper] - timebase[lower])
            return slope*(start_time - end_time)/2.0

#-------------------------------------------------------------------------------
#  Class to represent channels with a common time base.
#-------------------------------------------------------------------------------
class mpex_channel_group:
#-------------------------------------------------------------------------------
#  Initialize a mpex_channel_group.
#
#  Loads channels for the group.
#
#  param[inout] self       A mpex_channel_group instance.
#  param[in]    connection A connection to the MDSPlus server.
#  param[in]    group      Tree name of the group.
#  param[in]    channels   Channels in the group.
#-------------------------------------------------------------------------------
    def __init__(self, connection, group, channels):
        self.group = group
        self.channels = channels

        if 'timebase' in channels:
            self.timebase = mpex_channel(connection, '{}:{}'.format(self.group, self.channels['timebase']['name']), self.channels['timebase']['operations'])
            del channels['timebase']
        else:
            self.timebase = None

        if 'position' in channels:
            self.position = mpex_channel(connection, '{}:{}'.format(self.group, self.channels['position']['name']), self.channels['position']['operations'])
            del channels['position']
        else:
            self.position = None

        bad_keys = []
        for key in self.channels:
            try:
                self.channels[key] = mpex_channel(connection, '{}:{}'.format(self.group, key), self.channels[key]['operations'])
            except Exception as e:
                bad_keys.append(key)

        for key in bad_keys:
            print('No data for channel {}'.format(key))
            del self.channels[key]

#-------------------------------------------------------------------------------
#  Plot a mpex_channel_group:
#
#  Mpex channels can be 1, 2, or 3 dimensional + time.
#
#  param[in]    self     A mpex_channel_group instance.
#  param[in]    index    Index of the group to plot.
#  param[inout] axis     Area to plot data to.
#  param[in]    min_time Minimum plot time.
#  param[in]    max_time Maximum plot time.
#-------------------------------------------------------------------------------
    def plot(self, index, axis, min_time, max_time):
        key = list(self.channels)[index]
        self.channels[key].plot(axis, self.timebase, self.position, min_time, max_time)

#-------------------------------------------------------------------------------
#  Print a mpex_channel.
#
#  param[in] self A mpex_channel_group instance.
#-------------------------------------------------------------------------------
    def print(self):
        print(self.group)
        for key in self.channels:
            print('    ', end='')
            self.channels[key].print()

        if self.timebase != None:
            print('timebase : ', end='')
            self.timebase.print()
        if self.position != None:
            print('position : ', end='')
            self.position.print()

#-------------------------------------------------------------------------------
#  Get the number of channels
#
#  param[in] self A mpex_channel_group instance.
#  returns The number of channels in the group.
#-------------------------------------------------------------------------------
    def get_num_channels(self):
        return len(self.channels)

#-------------------------------------------------------------------------------
#  Get the min time.
#
#  param[in] self A mpex_channel_group instance.
#  returns The lowest time base.
#-------------------------------------------------------------------------------
    def get_min_time(self):
        if self.timebase != None:
            return numpy.ndarray.min(self.timebase.value)
        else:
            return None

#-------------------------------------------------------------------------------
#  Get the min time.
#
#  param[in] self A mpex_channel_group instance.
#  returns The lowest time base.
#-------------------------------------------------------------------------------
    def get_max_time(self):
        if self.timebase != None:
            return numpy.ndarray.max(self.timebase.value)
        else:
            return None

#-------------------------------------------------------------------------------
#  Slice a group.
#
#  param[in] self       A mpex_channel_group insrtance.
#  param[in] start_time Start time of the slice.
#  param[in] end_time   End time of the slice.
#-------------------------------------------------------------------------------
    def slice(self, start_time, end_time):
        group_slice = {}

        for key in self.channels:
            group_slice[key] = self.channels[key].slice(start_time, end_time, self.timebase.value)

        return group_slice

#-------------------------------------------------------------------------------
#  Get channel names.
#
#  param[in] self A mpex_channel_group instance.
#-------------------------------------------------------------------------------
    def get_channel_names(self):
        return self.channels

#-------------------------------------------------------------------------------
#  Class to represent a single shot.
#
#  Shots contain multiple channels.
#-------------------------------------------------------------------------------
class mpex_shot:
#-------------------------------------------------------------------------------
#  Initialize a mpex_shot.
#
#  This loads channels data for a shot.
#
#  param[inout] self       A mpex_shot instance.
#  param[in]    connection A connection to the MDSPlus server.
#  param[in]    number     Shot number for the data.
#  param[in]    config     Channel configureation.
#-------------------------------------------------------------------------------
    def __init__(self, connection, number, config):
        self.connection = connection
        self.number = number

        print('Opening mpex shot {}'.format(self.number))
        self.connection.get('TreeOpen("mpex",$)', self.number)

#  Present config from getting over written by forcing a copy.
        self.channels = copy.deepcopy(config)

        bad_keys = []
        for key in self.channels:
            try:
                self.channels[key] = mpex_channel_group(self.connection, key, self.channels[key])
            except Exception as e:
                bad_keys.append(key)

        for key in bad_keys:
            print('No timebase of position for group {}'.format(key))
            del self.channels[key]

#-------------------------------------------------------------------------------
#  Finialize mpex_shot before garbage collection.
#
#  param[inout] self A mpex_shot instance.
#-------------------------------------------------------------------------------
    def close(self):
        print('Closing mpex shot {}'.format(self.number))
        self.connection.get('TreeClose("mpex",$)', self.number)

#-------------------------------------------------------------------------------
#  Plot a mpex_shot:
#
#  Loop through each shot and present a syncronized graph of each channel.
#
#  param[in] self A mpex_shot instance.
#-------------------------------------------------------------------------------
    def plot(self):
        rows = self.get_num_channels()
        min_time = self.get_min_time()
        max_time = self.get_max_time()

        axes = matplotlib.pyplot.figure(figsize=(10, 13), constrained_layout=True).subplots(rows, 1)

        for index, axis in enumerate(axes):
            start = 0
            for key in self.channels:
                end = self.channels[key].get_num_channels()

                if start <= index and index < start + end:
                    self.channels[key].plot(index - start, axis, min_time, max_time)
                    break
                else:
                    start += end

        matplotlib.pyplot.show()

#-------------------------------------------------------------------------------
#  Print a mpex_channel.
#
#  param[in] self A mpex_shot instance.
#-------------------------------------------------------------------------------
    def print(self):
        for key in self.channels:
            self.channels[key].print()

#-------------------------------------------------------------------------------
#  Get the number of groups
#
#  param[in] self A mpex_shot instance.
#  returns The number of channels in the group.
#-------------------------------------------------------------------------------
    def get_num_channels(self):
        num = 0
        for key in self.channels:
            num += self.channels[key].get_num_channels()
        return num

#-------------------------------------------------------------------------------
#  Get the min time.
#
#  param[in] self A mpex_shot instance.
#  returns The lowest time base.
#-------------------------------------------------------------------------------
    def get_min_time(self):
        min_time = 1000000000
        for key in self.channels:
            group_min = self.channels[key].get_min_time()
            if group_min != None:
                min_time = min(min_time, group_min)
        return min_time

#-------------------------------------------------------------------------------
#  Get the min time.
#
#  param[in] self A mpex_shot instance.
#  returns The highest time base.
#-------------------------------------------------------------------------------
    def get_max_time(self):
        max_time = -1000000000
        for key in self.channels:
            group_max = self.channels[key].get_max_time()
            if group_max != None:
                max_time = max(max_time, self.channels[key].get_max_time())
        return max_time

#-------------------------------------------------------------------------------
#  Slice a shot
#
#  param[in] self       A mpex_shot instance.
#  param[in] start_time Start time of the slice.
#  param[in] end_time   End time of the slice.
#-------------------------------------------------------------------------------
    def slice(self, start_time, end_time):
        slice_data = {}
        for key in self.channels:
            slice_data[key] = self.channels[key].slice(start_time, end_time)

        return slice_data

#-------------------------------------------------------------------------------
#  Get channel names.
#
#  param[in] self A mpex_shot instance.
#-------------------------------------------------------------------------------
    def get_channel_names(self):
        names = []

        for key in self.channels:
            names += self.channels[key].get_channel_names()

        return names

#-------------------------------------------------------------------------------
#  Class to represent a time slice.
#-------------------------------------------------------------------------------
class mpex_time_slice:
#-------------------------------------------------------------------------------
#  Initialize a mpex_time_slice.
#
#  This loads channels data for a shot.
#
#  param[inout] self       A mpex_time_slice instance.
#  param[in]    shot_data  A connection to the MDSPlus server.
#  param[in]    start_time Start time of the slice.
#  param[in]    end_time   End time of the slice.
#-------------------------------------------------------------------------------
    def __init__(self, shot_data, start_time, end_time):
        self.shot_number = shot_data.number
        self.time = int(((end_time + start_time)/2.0)*1000)

        self.channels = shot_data.slice(start_time, end_time)

#-------------------------------------------------------------------------------
#  Print a mpex_channel.
#
#  param[in] self A mpex_time_slice instance.
#-------------------------------------------------------------------------------
    def print(self):
        for key in self.channels:
            print(key)
            print(self.channels[key])

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
                                     '--shotnumber',
                                     action='store',
                                     default=31527,
                                     type=int,
                                     dest='shotnumber',
                                     help='Shot number to grab mdsplus data from.',
                                     metavar='SHOTNUMBER')

    args = vars(command_line_parser.parse_args())

#  Remove empty arguments
    for key in [key for key in args if args[key] == None]:
        del args[key]

    with open(args['config'], 'r') as json_ref:
        config = json.load(json_ref)

    connection = MDSplus.Connection(args['mdsplusserver'])

    shot = mpex_shot(connection, args['shotnumber'], config)
    shot.print()

    slice = mpex_time_slice(shot, 4.0, 5.0)

    shot.plot()
    shot.close()
