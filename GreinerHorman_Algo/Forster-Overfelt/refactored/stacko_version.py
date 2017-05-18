# -*- coding: cp1252 -*-
'''
  The following is an adaptation of the above Greiner-Hormann* algorithm to deal
  with degenerate cases. The adaptation was briefly described by Liu et al.**  
  *Greiner, G. and Hormann K., Efficient Clipping of Arbitrary Polygons, ACM
  Trans. on Graphics, 17(2), 1998, pp.71-83
  **Liu, Y. K., Wang X. Q., Bao S. Z., Gombosi M., and Zalik B, An Algorithm for
  Polygon Clipping and for Determining Polygon Intersections and Unions, Comp. &
  Geo, 33, pp. 589-598, 2007
'''

# CLASSES

class Vertex(object):
    """Node in a circular doubly linked list.

    This class is almost exactly as described in the paper by Günther/Greiner.
    """

    def __init__(self, vertex, alpha=0.0, intersect=False, entry=True, checked=False, degen=False):
        if isinstance(vertex, Vertex):
            vertex = (vertex.x, vertex.y)
            # checked = True

        self.x, self.y = vertex     # point coordinates of the vertex
        self.next = None            # reference to the next vertex of the polygon
        self.prev = None            # reference to the previous vertex of the polygon
        self.neighbour = None       # reference to the corresponding intersection vertex in the other polygon
        self.entry = entry          # True if intersection is an entry point, False if exit
        self.degen = degen
        self.alpha = alpha          # intersection point's relative distance from previous vertex
        self.intersect = intersect  # True if vertex is an intersection
        self.checked = checked      # True if the vertex has been checked (last phase)

    def __repr__(self):
        if self.intersect:
            return "intersection%s"%str((self.x,self.y))
        else: return str((self.x,self.y))


class Polygon(object):
    """Manages a circular doubly linked list of Vertex objects that represents a polygon."""

    first = None

    def __iter__(self):
        """Iterator generator for this doubly linked list."""
        s = self.first
        while True:
            yield s
            s = s.__next__
            if s == self.first:
                return

    def add(self, vertex):
        """Add a vertex object to the polygon (vertex is added at the 'end' of the list")."""
        if not self.first:
            self.first = vertex
            self.first.next = vertex
            self.first.prev = vertex
        else:
            next = self.first
            prev = next.prev
            next.prev = vertex
            vertex.next = next
            vertex.prev = prev
            prev.next = vertex

    def replace(self, old, new):
        new.next = old.__next__
        new.prev = old.prev
        old.prev.next = new
        old.next.prev = new
        if old is self.first:
            self.first = new

    def insert(self, vertex, start, end):
        """Insert and sort a vertex between a specified pair of vertices.

        This function inserts a vertex (most likely an intersection point)
        between two other vertices (start and end). These other vertices
        cannot be intersections (that is, they must be actual vertices of
        the original polygon). If there are multiple intersection points
        between the two vertices, then the new vertex is inserted based on
        its alpha value.
        """
        curr = start
        while curr != end and curr.alpha < vertex.alpha:
            curr = curr.__next__

        vertex.next = curr
        prev = curr.prev
        vertex.prev = prev
        prev.next = vertex
        curr.prev = vertex

    @property
    def points(self):
        """Return the polygon's points as a list of tuples (ordered coordinates pair)."""
        p = []
        for v in self.iter():
            p.append((v.x, v.y))
        return p


# FUNCTIONS

def label_cases(current):
    """
    Deal with the cases:
    on/on, on/out, on/in, out/on, out/out, out/in, in/on, in/out, in/in
    """
    neighbour = current.neighbour
    #on/on
    if current.prev.intersect and current.next.intersect:
      #Determine what to do based on the neighbour
      #en tag is the opposite of the neighbour's en tag 
      if neighbour.prev.intersect and neighbour.next.intersect:
          current.intersect = False
          current.entry = True
          neighbour.entry = True
      elif neighbour.prev.intersect and not neighbour.next.entry:
          current.entry = False
      elif neighbour.prev.intersect and neighbour.next.entry:
          current.entry = True
      elif not neighbour.prev.entry and neighbour.next.intersect:
          current.entry = False
      elif not (neighbour.prev.entry or neighbour.next.entry):
          current.intersect = False
          current.entry = True
          neighbour.entry = False
      elif not neighbour.prev.entry and neighbour.next.entry:
          current.entry = False
      elif neighbour.prev.entry and neighbour.next.intersect:
          current.entry = True
      elif neighbour.prev.entry and not neighbour.next.entry:
          current.entry = True
      elif neighbour.prev.entry and neighbour.next.entry:
          current.intersect = False
          current.entry = False
          neighbour.entry = True
    #on/out
    elif current.prev.intersect and not current.next.entry:
        current.entry = False
    #on/in  
    elif current.prev.intersect and current.next.entry:
        current.entry = True
    #out/on  
    elif not current.prev.entry and current.next.intersect:
        current.entry = True
    #out/out  
    elif not (current.prev.entry or current.next.entry):
        if neighbour.prev.intersect and neighbour.next.intersect:
            current.intersect = False
            neighbour.entry = True
        elif neighbour.prev.entry == neighbour.next.entry:
            current.intersect = False
        else:
            if neighbour.prev.entry and not neighbour.next.entry:
                current.entry = True
            else:
                current.entry = False
    #out/in  
    elif not current.prev.entry and current.next.entry:
        current.entry = True
    #in/on
    elif current.prev.entry and current.next.intersect:
        current.entry = False
    #in/out
    elif current.prev.entry and not current.next.entry:
        current.entry = False
    #in/in
    elif current.prev.entry and current.next.entry:
        if neighbour.prev.intersect and neighbour.next.intersect:
            current.intersect = False
            neighbour.entry = False
        elif neighbour.prev.entry == neighbour.next.entry:
            current.intersect = False
        else:
            if neighbour.prev.entry and not neighbour.next.entry:
                current.entry = True
            else:
                current.entry = False

def testLocation(point, polygon):
    """
    Effective scanline test for the location of a point vis a vis a polygon.
    Returns either "in","on",or "out".
    Based on algorithm 7 from:
        Kai Horman and Alexander Agathos,
        "The point in polygon problem for arbitrary polygons".
        Computational Geometry: Theory and Applications, 
        Volume 20 Issue 3, November 2001
    """
    # begin
    if polygon.first.y == point.y and polygon.first.x == point.x:
        return "on" # vertex
    w = 0
    for v in polygon:
        if v.next.y == point.y:
            if v.next.x == point.x:
                return "on" # vertex
            else:
                if v.y == point.y and (v.next.x > point.x) == (v.x < point.x):
                    return "on" # edge
        # if crossing horizontal line
        if (v.y < point.y and v.next.y >= point.y)\
               or (v.y >= point.y and v.next.y < point.y):
            if v.x >= point.x:
                if v.next.x > point.x:
                    # modify w
                    if v.next.y > v.y: w += 1
                    else: w -= 1
                else:
                    det = (v.x - point.x) * (v.next.y - point.y) \
                        - (v.next.x - point.x) * (v.y - point.y)
                    if det == 0: return "on" # edge
                    # if right crossing
                    if (det > 0 and v.next.y > v.y)\
                       or (det < 0 and v.next.y < v.y):
                        # modify w
                        if v.next.y > v.y: w += 1
                        else: w -= 1
            else:
                if v.next.x > point.x:
                    det = (v.x - point.x) * (v.next.y - point.y) \
                        - (v.next.x - point.x) * (v.y - point.y)
                    if det == 0: return "on" # edge
                    # if right crossing
                    if (det > 0 and v.next.y > v.y)\
                       or (det < 0 and v.next.y < v.y):
                        # modify w
                        if v.next.y > v.y: w += 1
                        else: w -= 1
    if (w % 2) != 0:
        return "in"
    else:
        return "out"

def mark_locations(poly1, poly2):
    # looping both subj and clip, marking each vertex location
    for s in poly1:
        s.loc = testLocation(s, poly2)
    for c in poly2:
        c.loc = testLocation(c, poly1)

def insert_intersections(subj, clip):
    for s in subj: # for each vertex Si of subject polygon do
        if not s.intersect:
            for c in clip: # for each vertex Cj of clip polygon do
                if not c.intersect:
                    intersection = lines_intersect(s, s.__next__, c, c.__next__)
                    if intersection:
                        i, alphaS, alphaC = intersection
                        s_between = (0 < alphaS < 1)
                        c_between = (0 < alphaC < 1)
                        if s_between and c_between:
                            #both subj and clip intersect each other somewhere in the middle
                            iS = Vertex(i, alphaS, intersect=True, entry=False)
                            iC = Vertex(i, alphaC, intersect=True, entry=False)
                            subj.insert(iS, s, s.__next__)
                            clip.insert(iC, c, c.__next__)
                        else:
                            if s_between:
                                #subj line is touched by the start or stop point of a line from the clip polygon, so insert and mark that intersection as a degenerate
                                iS = Vertex(i, alphaS, intersect=True, entry=False, degen=True)
                                subj.insert(iS, s, s.__next__)
                            elif alphaS == 0:
                                #subj line starts at intersection, so mark the "degen"-flag, and replace vertex instead of inserting
                                iS = Vertex(i, alphaS, intersect=True, entry=False, degen=True)
                                subj.replace(s, iS)
                            elif alphaS == 1:
                                #subj line ends at intersection, so mark the "degen"-flag, and replace vertex instead of inserting
                                iS = Vertex(i, alphaS, intersect=True, entry=False, degen=True)
                                subj.replace(s.__next__, iS)
                            if c_between:
                                #clip line is touched by the start or stop point of a line from the subj polygon, so insert and mark that intersection as a degenerate
                                iC = Vertex(i, alphaC, intersect=True, entry=False, degen=True)
                                clip.insert(iC, c, c.__next__)
                            elif alphaC == 0:
                                #clip line starts at intersection, so mark the "degen"-flag, and replace vertex instead of inserting
                                iC = Vertex(i, alphaC, intersect=True, entry=False, degen=True)
                                clip.replace(c, iC)
                            elif alphaC == 1:
                                #clip line ends at intersection, so mark the "degen"-flag, and replace vertex instead of inserting
                                iC = Vertex(i, alphaC, intersect=True, entry=False, degen=True)
                                clip.replace(c.__next__, iC)
                            
                        iS.neighbour = iC
                        iC.neighbour = iS

#label intersections as entry or exit
def process_intersections(poly1, poly2):
    #cycle through first polygon and label intersections as en or ex
    current = poly1.first
    flag = True
    while flag:
        if current.intersect:
            pre = current.intersect
            label_cases(current)
            if current.interesct != pre: print(current)
            #Make sure current is still an intersection
            if current.intersect:
                label_cases(current.neighbour)
                #if the intersection is en/en or ex/ex
                if current.entry == current.neighbour.entry:
                    current.intersect = False
        if current == poly1.first:
            flag = False
        current = current.__next__ #move to the next point

def clip(subject, constraint):

    #make polygons
    subject = make_polygon(subject)
    constraint = make_polygon(constraint)
    clipped = []

    #prepping process
    mark_locations(subject, constraint) #label vertices as inside or outside
    insert_intersections(subject, constraint) #find intersections
    for s in subject: print(s)
    for c in constraint: print(c)
    process_intersections(subject, constraint) #label intersections and entry or exit and possibly remove

    flag = True #loop flag

    #set our current location to the first point in subject
    current = subject.first
    #loop through our polygon until we have found the first intersection
    while flag:
        current = current.__next__
        #Either an intersection has been found or no intersections found
        if current.intersect or current == subject.first:
            flag = False

    #reset our flag for the new loop
    flag = True
    #If a point lies outside of c and there was an intersection clip s
    if current.intersect:
        clipped.append((current.x,current.y))
        while flag:
            #Entry
            print("Hmm",current)
            if current.entry:
                current = current.__next__
                while not current.intersect:
                    clipped.append((current.x,current.y))
                    current = current.__next__
            #Exit
            else:
                current = current.prev
                while not current.intersect:
                    clipped.append((current.x,current.y))
                    current = current.prev

            #Check to see if we have completed a polygon
            if current.checked:
                #check if the polygon intersect at a point
                if len(clipped) is not 1:
                    #remove the last vertex because it is also the first 
                    del clipped[-1]
                #we have created our polygon so we can exit
                flag = False

            #change to the neighbour polygon since we are at a new intersection
            current.checked = True
            current = current.neighbour

    #Check if one polygon is contained wholly within the other
    elif contained(subject, constraint):
        clipped = subject.points
    elif contained(constraint, subject):
        clipped = constraint.points

    return clipped

def lines_intersect(s1, s2, c1, c2):
    """Same as intersect(), except returns
    intersection even if degenerate.
    """
    den = float( (c2.y - c1.y) * (s2.x - s1.x) - (c2.x - c1.x) * (s2.y - s1.y) )
    if not den:
        return None

    us = ((c2.x - c1.x) * (s1.y - c1.y) - (c2.y - c1.y) * (s1.x - c1.x)) / den
    uc = ((s2.x - s1.x) * (s1.y - c1.y) - (s2.y - s1.y) * (s1.x - c1.x)) / den

    if (0 <= us <= 1) and (0 <= uc <= 1):
        #subj and clip line intersect eachother somewhere in the middle
        #this includes the possibility of degenerates (edge intersections)
        x = s1.x + us * (s2.x - s1.x)
        y = s1.y + us * (s2.y - s1.y)
        return (x, y), us, uc
    else:
        return None

def contained(poly1, poly2):
    "test if poly1 is fully inside poly2"
    for p in poly1:
        if testLocation(p, poly2) == "out":
            return False
    return True

def make_polygon(points):
    poly = Polygon()
    for p in points:
        poly.add(Vertex(p))
    return poly

if __name__ == "__main__":
    """
    Test and visualize various polygon overlap scenarios.
    Visualization requires the pure-Python PyDraw library from
    https://github.com/karimbahgat/PyDraw
    """    
    import random
    subjpoly = [(0,0),(6,0),(6,6),(0,6),(0,0)]
    # normal intersections
    clippoly = [(4,4),(10,4),(10,10),(4,10),(4,4)] #simple overlap
    #clippoly = [(1,4),(3,8),(5,4),(5,10),(1,10),(1,4)] #jigzaw overlap
    #clippoly = [(7,7),(7,9),(9,9),(9,7),(7,7)] #smaller, outside
    #clippoly = [(2,2),(2,4),(4,4),(4,2),(2,2)] #smaller, inside
    #clippoly = [(-1,-1),(-1,7),(7,7),(7,-1),(-1,-1)] #larger, covering all
    #clippoly = [(-10,-10),(-10,-70),(-70,-70),(-70,-10),(-10,-10)] #larger, outside
    # degenerate intersections
    #clippoly = [(0,5),(6,4),(10,4),(10,10),(4,10),(0,5)] #degenerate, starts on edge intersection and goes inside
    #clippoly = [(5,6),(5.2,5.5),(5,5.4),(4.8,5.5)] #degenerate, starts on edge intersection and goes outside
    #clippoly = [(1,5),(6,4),(6,5),(10,4),(10,10),(4,10),(2,6),(1,6),(1,5)] #degenerate, hesitating to enter and exit
    #clippoly = [(1,5),(6,4),(6,5),(10,4),(10,10),(4,10),(2,6),(1.3,6),(1.6,6),(1,6),(1,5)] #degenerate, also multiple degens along shared line
    #clippoly = [(1,5),(6,4),(6,5),(10,4),(10,10),(4,10),(2,6),(1.5,5.7),(1,6),(0,6),(1,5)] #degenerate, back and forth on-out along shared line
    #clippoly = [(0,0),(6,0),(6,6),(0,6),(0,0)] #degenerate, perfect overlap
    #clippoly = [(1,0),(6,0),(6,6),(1,6),(1,0)] #degenerate, partial inside overlap
    #clippoly = [(0,6),(6,6),(6,10),(0,10),(0,6)] #degenerate, right next to eachother
    #clippoly = [(2,6),(6,6),(6,10),(2,10),(2,6)] #degenerate, partial right next to eachother
    # self intersecting polygons
    #clippoly = [(1,4),(3,8),(1,5),(5,4),(5,10),(1,10),(1,4)]
    # random polygon
    #clippoly = [(random.randrange(0,10),random.randrange(0,10)) for _ in xrange(10)] #random
    #run operation
    import time
    t = time.time()
    resultpoly = clip(subjpoly,clippoly)
    print("finished:",resultpoly,time.time()-t)
    import pydraw
    crs = pydraw.CoordinateSystem([-1,-1,11,11])
    img = pydraw.Image(400,400, crs=crs)
    img.drawpolygon(subjpoly, fillcolor=(222,0,0))
    img.drawpolygon(clippoly, fillcolor=(0,222,0))
    img.drawpolygon(resultpoly)
    img.drawgridticks(1,1)
    img.view()
