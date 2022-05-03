
from shapely.geometry import asShape, box as Box, mapping, LineString, MultiLineString,LinearRing
import fiona
import os
import numpy as np
import glob


def removeAloneSegments(lines):
    buffer = []
    [buffer.append(line) for line in lines if (line[1] - line[0]) > 2]
    return buffer

def queryCoast(coastLine, bbox):
    return ((shape, line) for shape, line in
            ((asShape(line["geometry"]), line) for line in coastLine)
            if shape.intersects(bbox))


def cutCoast(coastLine, bbox):
    for shape, line in coastLine:
        shape = shape.intersection(bbox)
        line["geometry"] = mapping(shape.intersection(bbox))
        yield line

class Converter:
    def __init__(self):
        self.open = []
        self.closedSizes = [0]
        self.openSizes = [0]
        self.closed = []

    def parse(self, data):
        self.open = []
        self.closed = []
        self.openSizes = [0]
        self.closedSizes = [0]
        for geo in data:
            if geo["type"] == "MultiLineString":
                for c in geo["coordinates"]:
                    self._push(c)
            else:
                self._push(geo["coordinates"])

    def _push(self, coords):
        if coords[0] == coords[-1]:
            self.closed.append(coords[:-1])
            self.closedSizes.append(self.closedSizes[-1] + len(coords) - 1)
        else:
            self.open.append(coords)
            self.openSizes.append(self.openSizes[-1] + len(coords))

    def toGeo(self, data):
        self.parse(data)

        coords = np.vstack(self.closed) if len(self.closed) > 0 else np.array([])

        sizes = self.closedSizes
        header = ["// closed curves"]
        indexes = ["IP = newp;", "IL = newl;", "IS = news;", "IF = newf;"]
        points = ["Point(IP + %d) = {%s, %s, 0}; //%d" % (i, c[0], c[1], i)
                  for i, c in enumerate(coords)]
        lines = ["Line(IL + %d) = {IP + %d : IP + %d, IP + %d};" % (
            i, sizes[i], sizes[i + 1] - 1, sizes[i])
                 for i in np.arange(len(sizes) - 1)] + ["//****lastClosed=IL+%d" % (len(sizes) - 1)]

        closedBuffer = header + indexes + points + lines

        coords = np.vstack(self.open)
        sizes = self.openSizes
        header = ["// open curves"]
        indexes = ["IP = newp;", "IL = newl;", "IS = news;", "IF = newf;"]
        points = ["Point(IP + %d) = {%s, %s, 0}; //%d" % (i, c[0], c[1], i)
                  for i, c in enumerate(coords)]
        lines = ["Line(IL + %d) = {IP + %d : IP + %d};" % (
            i, sizes[i], sizes[i + 1] - 1)
                 for i in np.arange(len(sizes) - 1)]
        openBuffer = header + indexes + points + lines

        return "\n".join(closedBuffer + openBuffer) + '\n'

    def toNPZ(self, data, outname):
        self.parse(data)
        np.savez(outname,closed=self.closed,closedSizes=self.closedSizes,open=self.open,openSizes=self.openSizes)
        print ('done')

    def toSHP(self, data, outname):
        self.parse(data)
        coords = np.vstack(self.closed) if len(self.closed) > 0 else np.array([])
        sizes = self.closedSizes
        lines = [(sizes[i], sizes[i + 1] ) for i in np.arange(len(sizes) - 1)]



        def checkSelfIntersection(poly):
            return poly.is_valid
        lines=removeAloneSegments(lines)
        polygons = []

        for line in lines:
            poly=[]
            for l in np.arange(line[0], line[1]):
                poly.append((coords[l][0], coords[l][1]))
            poly=LinearRing(poly)
            if checkSelfIntersection(poly):
                polygons.append(poly)
            else:
                poly = []
                for l in np.arange(line[0], line[1])[::-1]:
                    poly.append((coords[l][0], coords[l][1]))
                polygons.append(poly)


        coords = np.vstack(self.open)
        sizes = self.openSizes
        lines = [(sizes[i], sizes[i + 1]) for i in np.arange(len(sizes) - 1)]
        lines = removeAloneSegments(lines)
        # for line in lines:
        #     polygons.append(LinearRing([(coords[l][0],coords[l][1]) for l in np.arange(line[0],line[1])]))
        for line in lines:
            poly=[]
            for l in np.arange(line[0], line[1]):
                poly.append((coords[l][0], coords[l][1]))
            poly=LinearRing(poly)
            if checkSelfIntersection(poly):
                polygons.append(poly)
            else:
                poly = []
                for l in np.arange(line[1]-1, line[0],-1):
                    poly.append((coords[l][0], coords[l][1]))
                polygons.append(poly)

        print('converting on shapely Polygon')
        poly = MultiLineString(polygons)
        print('saving')
        # Define a polygon feature geometry with one attribute
        schema = {
            'geometry': 'MultiLineString',
            'properties': {'id': 'int'},
        }
        # Write a new Shapefile
        with fiona.open('%s.shp'%outname, 'w', 'ESRI Shapefile', schema) as c:
            c.write({
                'geometry': mapping(poly),
                'properties': {'id': 123},
            })
        print ('done')


class SegmentMerger:
    """
    Merge [Multi]Linestrings when terminations overlap. NO FORK ALOWED
    - closed strings are stored as they are presented
    - isolated strings are stored by terminations
    - if a string as a termination in common with a previously stored string, the two strings are merged together
    """

    def __init__(self):
        self.ends = {}
        self.outs = {}

    @staticmethod
    def join_lines(left, right):
        """ Join segments as 1D numpy array
        :param left: left-side segment
        :param right:  right-side segment
        :return: a new array containing the extended segment
        """
        return np.vstack([left[:-1], right[:]])

    def store(self, item):
        """ Store a segment (by head coordinate)
        :param item: array to store
        """
        self.outs[tuple(item[0])] = item

    def forget(self, item):
        """ Forget a previous stored item
        :param item: the item to forget about
        """
        del self.outs[tuple(item[0])]
        del self.ends[tuple(item[0])]
        del self.ends[tuple(item[-1])]

    def merge(self, geo):
        """ Merge/Store a new geometry
        :param geo: the segment(s) to merge or store
        :raise: TypeError if geo is a [Multi]LineString
        """
        if geo.type == "MultiLineString":
            [self.merge(i) for i in geo]
        elif geo.type == "LineString":
            coords = np.array(geo.coords)
            start = tuple(coords[0])
            end = tuple(coords[-1])
            if np.all(start == end):
                self.store(coords)
                return
            try:
                other = self.ends[start]
                self.forget(other)
                coords = self.join_lines(other, coords)
            except KeyError:
                pass
            try:
                other = self.ends[end]
                self.forget(other)
                coords = self.join_lines(coords, other)
            except KeyError:
                pass
            self.ends[tuple(coords[0])] = self.ends[tuple(coords[-1])] = coords
            self.store(coords)
        else:
            raise TypeError("only [Multi]LineString")


def asGeo(shapes, targetBox, simplify=None, updateCb=None):
    updateCb = updateCb if updateCb is not None else lambda *a: None
    features = list(cutCoast(queryCoast(shapes, targetBox), targetBox))
    updateCb("file read")
    shapes = [asShape(i["geometry"]) for i in features[:]]
    updateCb("shaped", len(shapes))
    merger = SegmentMerger()
    for i, shape in enumerate(shapes):
        merger.merge(shape)
        if i % (len(shapes) / 100.) < (i - 1) % (len(shapes) / 100.):
            updateCb("merging", int(i / (len(shapes) / 100.)))
    updateCb("merged")
    shapes = map(LineString, merger.outs.values())
    updateCb("joined")
    shapes = MultiLineString(shapes)
    updateCb("re-joined")
    if simplify is not None:
        shapes = shapes.simplify(simplify, True)
        updateCb("easier now")
    features = [mapping(shape) for shape in shapes]
    conv = Converter()
    return conv.toGeo(features)



def fromShp(box):

    box=Box(*box)

    mn = fiona.open(base)

    print ("open file")
    features = list(cutCoast(queryCoast(mn, box), box))

    # features = list(mn)
    print ("file read")
    shapes = [asShape(i["geometry"]) for i in features[:]]

    return shapes


def main():

    shapes=fromShp(box)

    print( "shaped", len(shapes))
    merger = SegmentMerger()
    for i, shape in enumerate(shapes):
        merger.merge(shape)
        if i % 1000 == 0:
            print (i / 1000),
    print ("-")
    shapes = list(map(LineString, merger.outs.values()))
    print ("joined")
    shapes = MultiLineString(shapes)

    print ("re-joined")

    shapes = shapes.simplify(simply, True)
    #print "easier now"
    features = [mapping(shape) for shape in shapes]

    conv = Converter()

    if outputType=='geo':
        print ("2geo")
        out = conv.toGeo(features)
        open(f'{outName}_simpl{simply}.geo', "w").write(out)
    elif outputType=='shp':
        conv.toSHP(features, f'{outName}_simpl{simply}')
    elif outputType=='npz':
        conv.toNPZ(features, f'{outName}.npz')
    else:
        exit(f'{outputType} not available')







outName = 'Sicily_N_CMED'
#base = "/Users/scausio/Dropbox (CMCC)/data/coastline/coastline/coastlinesOpenstreetMap/lines.shp"
#base='/Users/scausio/Dropbox (CMCC)/data/shapefile/gshhg-shp-2.3.6/GSHHS_shp/f/GSHHS_f_L1.shp'
base='/Users/scausio/Dropbox (CMCC)/data/coastline/coastline/gshhsCoastlineGlobaleHD/gshhsLines.shp'
box =  [4.69,28.47, 24.85,47.10]    # lon_min, lat_min,lon_max, lat_max

simply=0
outputType='shp' # available types are ['shp','npz', 'geo']

if __name__ == "__main__":
    main()
