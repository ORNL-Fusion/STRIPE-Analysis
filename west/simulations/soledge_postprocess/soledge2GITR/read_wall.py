from shapely.geometry import Point, Polygon

class Surface:
    def __init__(self, filename, geom_type):
        # read surface file
        with open(filename, 'r') as f:
            # skip first two lines
            data = f.readlines()[2:]
            # read the number of points and lines
            num_points = int(data[0].split()[0])
            num_lines = int(data[1].split()[0])
            print('Number of points: ', num_points)
            print('Number of lines: ', num_lines)
            # read points
            self.points = {}
            for line in data[5:num_points]:
                # print(line)
                point_id, *point_coords = map(float, line.split())
                self.points[int(point_id)] = tuple(point_coords)
            # read lines
            self.lines = []
            # self.material = {}
            for line in data[num_points+8:]:
                line_parts = line.split()
                materials= line_parts[-1]
                line_id= int(line_parts[0])
                line_points = [int(line_parts[1]), int(line_parts[2])]
                self.lines.append((int(line_id), str(materials), line_points))
        # create Polygon object from points
        points_list = [list(p) for p in self.points.values()]
        self.polygon = Polygon(points_list)
        # set surface ID
        self.id = filename.split('.')[0]
        
        
