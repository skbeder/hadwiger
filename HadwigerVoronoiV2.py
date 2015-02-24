import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
from scipy.spatial import distance as myDist
import random

class HadwigerVoronoi:

    def voronoi_finite_polygons_2d(self, vor, radius=None):
        """
        Reconstruct infinite voronoi regions in a 2D diagram to finite
        regions.

        Parameters
        ----------
        vor : Voronoi
            Input diagram
        radius : float, optional
            Distance to 'points at infinity'.

        Returns
        -------
        regions : list of tuples
            Indices of vertices in each revised Voronoi regions.
        vertices : list of tuples
            Coordinates for revised Voronoi vertices. Same as coordinates
            of input vertices, with 'points at infinity' appended to the
            end.

        """

        if vor.points.shape[1] != 2:
            raise ValueError("Requires 2D input")

        new_regions = []
        new_vertices = vor.vertices.tolist()

        center = vor.points.mean(axis=0)
        if radius is None:
            radius = vor.points.ptp().max()

        # Construct a map containing all ridges for a given point
        all_ridges = {}
        for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
            all_ridges.setdefault(p1, []).append((p2, v1, v2))
            all_ridges.setdefault(p2, []).append((p1, v1, v2))

        # Reconstruct infinite regions
        for p1, region in enumerate(vor.point_region):
            vertices = vor.regions[region]

            if all(v >= 0 for v in vertices):
                # finite region
                new_regions.append(vertices)
                continue

            #try:
            # reconstruct a non-finite region
            ridges = all_ridges[p1]
            new_region = [v for v in vertices if v >= 0]

            for p2, v1, v2 in ridges:
                if v2 < 0:
                    v1, v2 = v2, v1
                if v1 >= 0:
                    # finite ridge: already in the region
                    continue

                # Compute the missing endpoint of an infinite ridge

                t = vor.points[p2] - vor.points[p1] # tangent
                t /= np.linalg.norm(t)
                n = np.array([-t[1], t[0]])  # normal

                midpoint = vor.points[[p1, p2]].mean(axis=0)
                direction = np.sign(np.dot(midpoint - center, n)) * n
                far_point = vor.vertices[v2] + direction * radius

                new_region.append(len(new_vertices))
                new_vertices.append(far_point.tolist())

            # sort region counterclockwise
            vs = np.asarray([new_vertices[v] for v in new_region])
            c = vs.mean(axis=0)
            angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
            new_region = np.array(new_region)[np.argsort(angles)]

            # finish
            new_regions.append(new_region.tolist())
            #except:
            #    print "MAJOR UNKNOWN ERROR :-("

        return new_regions, np.asarray(new_vertices)

    def colorize(self):
        # make up data points
        #np.random.seed(1234)


        # compute Voronoi tesselation

        vor = Voronoi(self.points)

        # plot
        regions, vertices = self.voronoi_finite_polygons_2d(vor, 1)

        self.regions = regions
        self.vertices = vertices




    def showGraph(self):
         # colorize
        colors = ['b','g','r','c','y','k', '0.75']

        # for i in xrange(self.numcolors):
        #     for j in xrange(self.numinstance):
        for pointIndex in xrange(len(self.points)):
            point = self.points[pointIndex]
            color = self.reversecolormap[pointIndex]
            polygon = self.vertices[self.regions[pointIndex]]
            plt.fill(*zip(*polygon), facecolor = colors[color], alpha=0.8)      

        '''
        for region in regions:
            polygon = vertices[region]
            plt.fill(*zip(*polygon), facecolor = colors[random.randint(0,5)], alpha=0.8)
        '''

        for i in xrange(len(self.points)):
            plt.plot(self.points[i][0],self.points[i][1], 'kx')
        #plt.plot(points, 'kx')
        #plt.plot(points[:,0], points[:,1], 'kx')
        
        #plt.xlim(vor.min_bound[0] - 0.1, vor.max_bound[0] + 0.1)
        #plt.ylim(vor.min_bound[1] - 0.1, vor.max_bound[1] + 0.1)
        plt.xlim(-1*self.wallDistance,self.wallDistance)
        plt.ylim(-1*self.wallDistance,self.wallDistance)

        plt.show()



    def getInternalDistance(self,index):


        maxDist = 0.0
        for point1 in self.regions[index]:
            for point2 in self.regions[index]:
                distance = myDist.euclidean(self.vertices[point1],self.vertices[point2])
                if distance > maxDist:
                    maxDist = distance



        return maxDist


    def getAdjSet(self, index):
        currRegion = self.regions[index]
        adjSet = []
        for i in xrange(len(self.regions)):
            for point in self.regions[i]:
                if point in currRegion and i != index:
                    adjSet.append(i)

        return adjSet

    def moveAdj(self,adjSet,index):
        for point in adjSet:
            if self.reversecolormap[index] != self.reversecolormap[point]:
                if myDist.euclidean(self.vertices[point],self.vertices[index]) > self.internalThresh*self.internalMax:


                    self.points[point][0] = self.internalMoveFactor*self.points[point][0] + (1-self.internalMoveFactor)*self.points[index][0] 
                    self.points[point][1] = self.internalMoveFactor*self.points[point][1] + (1-self.internalMoveFactor)*self.points[index][1] 




    def internalMotion(self, index):

        internalDistance = self.getInternalDistance(index)


        if internalDistance > self.internalMax:
            adjSet = self.getAdjSet(index)
            self.moveAdj(adjSet, index)



    def pointsTooClose(self, point, index):
        minDistance = 1000000
        for vertex1 in self.regions[point]:
            for vertex2 in self.regions[index]:
                distance = myDist.euclidean(self.vertices[vertex1],self.vertices[vertex2])
                if distance < minDistance:
                    minDistance = distance
        return minDistance

    def forceWithinWalls(self, index):
        if self.points[index][0]>self.wallDistance:
            self.points[index][0] = self.wallDistance
        if self.points[index][0] < -1*self.wallDistance:
            self.points[index][0] = -1*self.wallDistance 
        if self.points[index][1]>self.wallDistance:
            self.points[index][1] = self.wallDistance
        if self.points[index][1] < -1*self.wallDistance:
            self.points[index][1] = -1*self.wallDistance      

    def movePointExternal(self,point, index):
        self.points[index][0] = (1+self.externalMoveFactor)*self.points[index][0] + (-1*self.externalMoveFactor)*self.points[point][0] 
        self.points[index][1] = (1+self.externalMoveFactor)*self.points[index][1] + (-1*self.externalMoveFactor)*self.points[point][1] 
        self.forceWithinWalls(index)


    def externalMotion(self,index):
        color = self.reversecolormap[index]
        colorPoints = self.colormap[color]

        for point in colorPoints:
            if point != index:
                distance = self.pointsTooClose(point, index)
                if(distance == 0):
                    if(self.numIterationsPassed > self.numIterationsBeforeInsert):
                        self.addShapeBetweenAdjacentColors(point,index)
                elif (distance < self.externalMin):
                    if myDist.euclidean(self.vertices[point],self.vertices[index]) < self.externalThresh*self.externalMin:
                        self.movePointExternal(point,index)


    def findBestColor(self, point, invalidColor):
        maxDistance = -1
        colorToReturn = -1
        for i in xrange(self.numcolors):
            if i != invalidColor:
                minColorDist = 1000000
                for colorPoint in self.colormap[i]:
                    for vertex in self.regions[colorPoint]:
                        dist = myDist.euclidean(self.vertices[vertex],point)
                        if dist < minColorDist:
                            minColorDist = dist
                if minColorDist > maxDistance:
                    colorToReturn = i
                    maxDistance = minColorDist
        return colorToReturn

    def addShapeBetweenAdjacentColors(self, point1, point2):
        xcoord = (self.points[point1][0]+self.points[point2][0])/2
        ycoord = (self.points[point1][1]+self.points[point2][1])/2
        newpoint = []
        newpoint.append(xcoord)
        newpoint.append(ycoord)

        for point in self.points:
            if newpoint[0] == point[0] and newpoint[1] == point[1]:
                print "IN HERE"
                return

        color = self.reversecolormap[point1]

        newcolor = self.findBestColor(newpoint, color)

        self.points.append(newpoint)
        self.reversecolormap[len(self.points)-1] = newcolor
        self.colormap[newcolor].append(len(self.points)-1)

        self.regions.append([len(self.regions)])
        np.append(self.vertices,newpoint)


    def magnetMotion(self):
        #for colorset in self.points:
        #    for point in colorset:
        #        point[0] += random.random()/20
        for i in xrange(len(self.points)):
            self.internalMotion(i)
            self.externalMotion(i)




    def setup(self):
        self.numIterationsPassed = 0
        self.numIterationsBeforeInsert = 10
        self.internalMax = 1
        self.externalMin = 1

        #self.internalMin = 0.01

        self.numcolors = 7
        self.numinstance = 17

        self.internalMoveFactor = 0.95
        self.externalMoveFactor = 0.05

        self.internalThresh = 0.1
        self.externalThresh = 4

        self.wallDistance = 10

        colormap = {}
        reversecolormap = {}


        points = []
        for i in xrange(self.numcolors):
            currmap = []
            for j in xrange(self.numinstance):
                point = []
                point.append(random.random()*self.wallDistance)
                point.append(random.random()*self.wallDistance) 
                #point.append(float((1+i*self.numinstance+j))/18 + random.random()/20)
                #point.append(float((1+i*self.numinstance+j))/18 + random.random()/20)
                points.append(point)
                #colorset.append((random.random(),random.random()))


                currmap.append((i*self.numinstance+j))
                reversecolormap[i*self.numinstance+j] = i


            colormap[i] = currmap
        
        self.colormap = colormap
        self.reversecolormap = reversecolormap

        self.points = points




    def hadwiger(self):
        self.setup()
        self.colorize()
        for i in xrange(600):
            print "Numpoints: ", len(self.points)
            print "iteration ", i
            self.magnetMotion()
            self.colorize()
            if i%5==0:
                self.showGraph()
            self.numIterationsPassed += 1
        self.colorize()
        self.showGraph()

tester = HadwigerVoronoi()
tester.hadwiger()
#HadwigerVoronoi(setup())