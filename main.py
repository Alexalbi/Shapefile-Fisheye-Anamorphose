# pip install pyshp
# pip install matplotlib
import inspect
import shapefile
from math import*
import matplotlib.pyplot as plt

# sets the x and y constants for readability
global x, y
x, y = 0, 1

# test file location
testpath = "tests/"
testfilename = "grid.shp"

def switchList(list):
    # switches the formatting of a list of points
    xList, yList = [list[i][x] for i in range(len(list))], [list[i][y] for i in range(len(list))]

    return xList, yList


def plotShapes(shapes, points=True):
    # plots a single shape
    for shape in shapes:
        if points:
            list = switchList(shape.points)
        else:
            list = switchList(shape)
        plt.plot(list[0], list[1], 'black')
    plt.show()


def changeTransformScale(function, newAtOrigin, newAtBound):
    # scales a manotously decreasing function on [0, 1] to vary between newAtOrigin and newAtBound
    oldAtOrigin, oldAtBound = function(0), function(1)
    ratio = (newAtOrigin - newAtBound) / (oldAtOrigin - oldAtBound)
    return lambda x: (function(x) - oldAtBound) * ratio + newAtBound


def plotTransform(function, ax1, n=100):
    t = [i/(n - 1) for i in range(n)]
    deform = [function(ti) for ti in t]
    tdeform = [t[i] * deform[i] for i in range(n)]
    ax1.plot(t, deform, color='blue', label='transform function')
    ax2 = ax1.twinx()
    ax2.plot(t, tdeform, color='red', label='applied transformation')
    ax1.legend()
    ax2.legend()


def getVector(pointA, pointB):
    # returns the vector pointing from A to B
    return [(pointB[x] - pointA[x]), (pointB[y] - pointA[y])]


def computeDistance(pointA, pointB, returnVector=False):
    # returns the distance between A and B and the associated vector
    vector = getVector(pointA, pointB)
    distance = sqrt(vector[x] ** 2 + vector[y] ** 2)
    if returnVector:
        return distance, vector
    else:
        return distance


def transformShape(origin, maxDistance, points, function):
    newPoints = []

    for point in points:
        # gets the current point distance to origin
        distance, vector = computeDistance(origin, point, True)

        # if the transform is anisotropic
        if type(maxDistance)==list:
            normalized = [0, 0]
            # normalizes distance depending on the 4 directions
            if vector[x] > 0:
                normalized[x] = distance / maxDistance[0]
            else:
                normalized[x] = distance / maxDistance[1]
            if vector[y] > 0:
                normalized[y] = distance / maxDistance[2]
            else:
                normalized[y] = distance / maxDistance[3]
            # computes the transform value for the point
            transform = [function(normalized[x]), function(normalized[y])]
            # scales the vector
            vector = [transform[x] * vector[x], transform[y] * vector[y]]
        else:
            # normalizes distance
            normalized = distance / maxDistance
            if normalized > 1:
                normalized = 1
            # computes the transform value for the point
            transform = function(normalized)
            # scales the vector
            vector = [transform * vector[x], transform * vector[y]]

        # computes the new point position
        newPoint = [origin[x] + vector[x], origin[y] + vector[y]]
        newPoints.append(newPoint)

    return newPoints


def computeTransform(origin, function, maxDistance, shapeFileR, shapeFileW):
    if (type(maxDistance) == list):
        if type(origin[0] == list):
            assert len(maxDistance) == len(origin)
        else:
            assert len(maxDistance) == 4
            origin = [origin]
    else:
        if type(origin[0] != list):
            maxDistance = [maxDistance]
            origin = [origin]

    for shapeRec in shapeFileR.iterShapeRecords():
        # computes the transform of the shape
        newShape = shapeRec.shape.points
        for i in range(len(origin)):
            newShape = transformShape(origin[i], maxDistance[i], newShape, function)

        # copies the shape's record
        shapeFileW.record(*shapeRec.record)

        # uses a unique write method for every shape type described
        if shapeRec.shape.shapeType == shapefile.POINT:
            shapeFileW.point(newShape[0][x],newShape[0][y])
        elif shapeRec.shape.shapeType == shapefile.POLYLINE:
            shapeFileW.line([newShape])
        elif shapeRec.shape.shapeType == shapefile.POLYGON:
            shapeFileW.poly([newShape])
        elif shapeRec.shape.shapeType == shapefile.MULTIPOINT:
            shapeFileW.multipoint([newShape])
        elif shapeRec.shape.shapeType == shapefile.POINTZ:
            shapeFileW.pointz([newShape])
        elif shapeRec.shape.shapeType == shapefile.POLYLINEZ:
            shapeFileW.linez([newShape])
        elif shapeRec.shape.shapeType == shapefile.POLYGONZ:
            shapeFileW.polyz([newShape])
        elif shapeRec.shape.shapeType == shapefile.MULTIPOINTZ:
            shapeFileW.multipointz([newShape])
        elif shapeRec.shape.shapeType == shapefile.POINTM:
            shapeFileW.pointm([newShape])
        elif shapeRec.shape.shapeType == shapefile.POLYLINEM:
            shapeFileW.linem([newShape])
        elif shapeRec.shape.shapeType == shapefile.POLYGONM:
            shapeFileW.polym([newShape])
        elif shapeRec.shape.shapeType == shapefile.MULTIPOINTM:
            shapeFileW.multipointm([newShape])
        elif shapeRec.shape.shapeType == shapefile.MULTIPATCH:
            shapeFileW.multipatch([newShape])


def computeFisheye(origin, filenames, inpath, outpath, transform, maxDistance, display=False):
    if display:
        fig, ax1 = plt.subplots()

        fstring = str(inspect.getsourcelines(transform)[0])
        fstring = fstring.strip("['\\n']").split(" = ")[1].split(" : ")[1]
        plt.title("f(x) = " + str(fstring) + " | zoom from " + str(transformOrigin) + " to " + str(transformBound))

        transform = changeTransformScale(transform, transformOrigin, transformBound)
        plotTransform(transform, ax1)
        fig.tight_layout()
        plt.show()

    inPaths = [inpath + filename for filename in filenames]
    outPaths = [outpath + filename for filename in filenames]

    for i in range(len(inPaths)):
        # create shapefile reader
        shapeFileR = shapefile.Reader(inPaths[i])
        # create shapefile writer using same type as reader
        shapeFileType = shapeFileR.shapeType
        shapeFileW = shapefile.Writer(outPaths[i], shapeType=shapeFileType)
        shapeFileW.autoBalance = 1
        # copy the existing fields
        shapeFileW.fields = shapeFileR.fields[1:]

        # compute shapes transform
        computeTransform(origin, transform, maxDistance, shapeFileR, shapeFileW)

        # save the new shapefile
        shapeFileW.close()


def displayShapefile(filename):
    shapeFileR = shapefile.Reader(filename)
    shapes = shapeFileR.shapes()
    plotShapes(shapes)


def testFisheye(origin, transform, maxDistance):
    computeFisheye(origin, [testfilename], testpath, testpath+"transformed_", transform, maxDistance, True)
    displayShapefile(testpath+"transformed_"+testfilename)

# center coordinates of the "fisheye" transform
origin = [[842440, 6519200], [844710, 6517720]]

# transform function, assumes it is defined and monotonously decreasing on [0; 1]
# note: x*transform(x) should also be monotonous on [0; 1]
transform = lambda x: -log(0.2*x+0.1)

# transform factors/zoom at origin and bound
transformOrigin = 1.5
transformBound = 1

# radius of the fisheye effect
maxDistance = [2500.0, 1200]

# shapefile path - .shp/.dbf/... library does not care about file extensions
#path = "base_shapes_Lyon/"
#filenames = ["bati_aroundLyon.shp", "bridges.shp", "cimetiere_aroundLyon.shp", "hydro_Lyon.shp", "metro_lines.shp",
#             "metro_stations.shp", "parcs_aroundLyon.shp", "places_aroundLyon.shp", "railway_lines_aroundLyon.shp", "streets_name.shp",
#           "terrain_sport_aroundLyon.shp", "tram_lines.shp", "tram_stations.shp", "veget_aroundLyon.shp"]


#computeFisheye(origin, filenames, path, path+"transformed_", transform, maxDistance)
testFisheye(origin, transform, maxDistance)


# todo refactor for clarity
# todo optimize for fisheye bound based chunks
# todo cleanup method for use outside of file