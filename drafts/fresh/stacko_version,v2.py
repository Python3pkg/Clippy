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
            s = s.next
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
        new.next = old.next
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
            curr = curr.next

        vertex.next = curr
        prev = curr.prev
        vertex.prev = prev
        prev.next = vertex
        curr.prev = vertex

    @property
    def points(self):
        """Return the polygon's points as a list of tuples (ordered coordinates pair)."""
        p = []
        for v in self:
            p.append((v.x, v.y))
        return p


# FUNCTIONS

def label_cases(current, mode):
    """
    Deal with the cases:
    on/on, on/out, on/in, out/on, out/out, out/in, in/on, in/out, in/in
    """
    #prep
    if mode == "intersect":
        s_entry = True
        c_entry = True
    elif mode == "difference":
        s_entry = False
        c_entry = True
    elif mode == "union":
        s_entry = False
        c_entry = False
    neighbour = current.neighbour
    #on/on
    if current.prev.intersect and current.next.intersect:
      #Determine what to do based on the neighbour
      #en tag is the opposite of the neighbour's en tag 
      if neighbour.prev.intersect and neighbour.next.intersect:
          current.intersect = False
          current.entry = s_entry
          neighbour.entry = c_entry
      elif neighbour.prev.intersect and not neighbour.next.entry:
          current.entry = not s_entry
      elif neighbour.prev.intersect and neighbour.next.entry:
          current.entry = s_entry
      elif not neighbour.prev.entry and neighbour.next.intersect:
          current.entry = not s_entry
      elif not (neighbour.prev.entry or neighbour.next.entry):
          current.intersect = False
          current.entry = s_entry
          neighbour.entry = not c_entry
      elif not neighbour.prev.entry and neighbour.next.entry:
          current.entry = not s_entry
      elif neighbour.prev.entry and neighbour.next.intersect:
          current.entry = s_entry
      elif neighbour.prev.entry and not neighbour.next.entry:
          current.entry = s_entry
      elif neighbour.prev.entry and neighbour.next.entry:
          current.intersect = False
          current.entry = not s_entry
          neighbour.entry = c_entry
    #on/out
    elif current.prev.intersect and not current.next.entry:
        current.entry = not s_entry
    #on/in  
    elif current.prev.intersect and current.next.entry:
        current.entry = s_entry
    #out/on  
    elif not current.prev.entry and current.next.intersect:
        current.entry = s_entry
    #out/out  
    elif not (current.prev.entry or current.next.entry):
        if neighbour.prev.intersect and neighbour.next.intersect:
            current.intersect = False
            neighbour.entry = c_entry
        elif neighbour.prev.entry == neighbour.next.entry:
            current.intersect = False
        else:
            if neighbour.prev.entry and not neighbour.next.entry:
                current.entry = s_entry
            else:
                current.entry = not s_entry
    #out/in  
    elif not current.prev.entry and current.next.entry:
        current.entry = s_entry
    #in/on
    elif current.prev.entry and current.next.intersect:
        current.entry = not s_entry
    #in/out
    elif current.prev.entry and not current.next.entry:
        current.entry = not s_entry
    #in/in
    elif current.prev.entry and current.next.entry:
        if neighbour.prev.intersect and neighbour.next.intersect:
            current.intersect = False
            neighbour.entry = not c_entry
        elif neighbour.prev.entry == neighbour.next.entry:
            current.intersect = False
        else:
            if neighbour.prev.entry and not neighbour.next.entry:
                current.entry = s_entry
            else:
                current.entry = not s_entry

def test_location(point, polygon):
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
        s.loc = test_location(s, poly2)
    for c in poly2:
        c.loc = test_location(c, poly1)

def insert_intersections(subj, clip):
    s_intsecs = []
    c_intsecs = []
    for s in subj: # for each vertex Si of subject polygon do
        for c in clip: # for each vertex Cj of clip polygon do
            intersection = lines_intersect(s, s.next, c, c.next)
            if intersection:
                i, alphaS, alphaC = intersection

                iS = Vertex(i, alphaS, intersect=True, entry=False)
                iC = Vertex(i, alphaC, intersect=True, entry=False)

                iS.neighbour = iC
                iC.neighbour = iS
                
                s_intsecs.append( (iS, s, s.next) )
                c_intsecs.append( (iC, c, c.next) )
                        
        for iS,s,s_next in s_intsecs:
            subj.insert(iS, s, s_next)
        for iC,c,c_next in c_intsecs:
            clip.insert(iC, c, c_next)

        # collect into lists and relink to fix any insertion link issues
        slist = [s for s in subj]
        subj.first = None
        for s in slist:
            subj.add(s)
        clist = [c for c in clip]
        clip.first = None
        for c in clist:
            clip.add(c)

#label intersections as entry or exit
def process_intersections(poly1, poly2, mode):
    #cycle through first polygon and label intersections as en or ex
    current = poly1.first
    flag = True
    while flag:
        if current.intersect:
            pre = current.intersect
            label_cases(current, mode)
            if current.intersect != pre: print current
            #Make sure current is still an intersection
            if current.intersect:
                label_cases(current.neighbour, mode)
                #if the intersection is en/en or ex/ex
                if current.entry == current.neighbour.entry:
                    current.intersect = False
        if current == poly1.first:
            flag = False
        current = current.next #move to the next point

def clip(subject, constraint, mode):

    #make polygons
    subject = make_polygon(subject)
    constraint = make_polygon(constraint)
    clipped = []

    #prepping process
    mark_locations(subject, constraint) #label vertices as inside or outside
    insert_intersections(subject, constraint) #find intersections
    for s in subject: print s
    for c in constraint: print c
    process_intersections(subject, constraint, mode) #label intersections and entry or exit and possibly remove

    flag = True #loop flag

    #set our current location to the first point in subject
    current = subject.first
    #loop through our polygon until we have found the first intersection
    while flag:
        current = current.next
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
            print "Hmm",current
            if current.entry:
                current = current.next
                while not current.intersect:
                    clipped.append((current.x,current.y))
                    current = current.next
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

    clipped = [(clipped,[])] #naah, forces multipoly

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
        if test_location(p, poly2) == "out":
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
    
    subjpoly = [(0,0),(6,0),(6,6),(0,6),(0,0)]
    
    # normal intersections
    testpolys_normal = {"simple overlap":
                        [(4,4),(10,4),(10,10),(4,10),(4,4)],
                        "jigzaw overlap":
                        [(1,4),(3,8),(5,4),(5,10),(1,10),(1,4)],
                        "smaller, outside":
                        [(7,7),(7,9),(9,9),(9,7),(7,7)],
                        "smaller, inside":
                        [(2,2),(2,4),(4,4),(4,2),(2,2)],
                        "larger, covering all":
                        [(-1,-1),(-1,7),(7,7),(7,-1),(-1,-1)],
                        "larger, outside":
                        [(-10,-10),(-10,-70),(-70,-70),(-70,-10),(-10,-10)]
                        }

    # degenerate intersections
    testpolys_degens = {"degenerate, starts on edge intersection and goes inside":
                        [(0,5),(6,4),(10,4),(10,10),(4,10),(0,5)],
                        "degenerate, starts on edge intersection and goes outside":
                        [(5,6),(5.2,5.5),(5,5.4),(4.8,5.5)],
                        "degenerate, hesitating to enter and exit":
                        [(1,5),(6,4),(6,5),(10,4),(10,10),(4,10),(2,6),(1,6),(1,5)],
                        "degenerate, also multiple degens along shared line":
                        [(1,5),(6,4),(6,5),(10,4),(10,10),(4,10),(2,6),(1.3,6),(1.6,6),(1,6),(1,5)],
                        "degenerate, back and forth on-out along shared line":
                        [(1,5),(6,4),(6,5),(10,4),(10,10),(4,10),(2,6),(1.5,5.7),(1,6),(0,6),(1,5)]
                        }
    
    # nextto/almost copy special cases
    testpolys_nextto_almostsame = {"degenerate, perfect overlap":
                                   [(0,0),(6,0),(6,6),(0,6),(0,0)],
                                   "degenerate, partial inside overlap":
                                   [(1,0),(6,0),(6,6),(1,6),(1,0)],
                                   "degenerate, right next to eachother":
                                   [(0,6),(6,6),(6,10),(0,10),(0,6)],
                                   "degenerate, partial right next to eachother":
                                   [(2,6),(6,6),(6,10),(2,10),(2,6)]
                                   }
    
    #run operation

    import os
    import time
    import pydraw

##    # test geo
##    
##    def test_draw(testname, subjpoly, clippoly, mode):
##        t = time.time()
##        #print testname, mode
##        resultpolys = clip_polygon(subjpoly,clippoly,mode)
##        print "finished:",len(resultpolys),time.time()-t
##        print "start",str(resultpolys)[:100]
##        print "end",str(resultpolys)[-100:]
##        crs = pydraw.CoordinateSystem([0,80,45,50])
##        img = pydraw.Image(300,300, crs=crs)
##        img.drawpolygon(subjpoly, fillcolor=(222,0,0,111))
##        img.drawpolygon(clippoly, fillcolor=(0,222,0,111))
##        for ext,holes in resultpolys:
##            img.drawpolygon(ext,holes)
##        img.drawgridticks(10,10)
##        img.save("test_output/"+testname+"-"+mode+".png")
##
##    import pygeoj
##    world = pygeoj.load("cshapes.geo.json")
##    norw = next(cntr.geometry.coordinates[0][0] for cntr in world if cntr.properties["CNTRY_NAME"] == "Norway")
##    swed = next(cntr.geometry.coordinates[0][0] for cntr in world if cntr.properties["CNTRY_NAME"] == "Sweden")
##    test_draw("norway-sweden", norw, swed, "intersect")
##    
##    breakonpurpose

    # test basics

    def test_draw(testname, subjpoly, clippoly, mode):
        t = time.time()
        #print testname, mode
        resultpolys = clip(subjpoly,clippoly,mode)
        #print "finished:",resultpolys,time.time()-t
        crs = pydraw.CoordinateSystem([-1,-1,11,11])
        img = pydraw.Image(300,300, crs=crs)
        img.drawpolygon(subjpoly, fillcolor=(222,0,0))
        img.drawpolygon(clippoly, fillcolor=(0,222,0))
        for ext,holes in resultpolys:
            img.drawpolygon(ext,holes)
        img.drawgridticks(1,1)
        img.save("test_output/"+testname+"-"+mode+".png")

    if not os.path.lexists("test_output"): os.mkdir("test_output")
    for testname,testclip in testpolys_normal.items():
        for mode in ("intersect","union","difference"):
            test_draw(testname, subjpoly, testclip, mode)
    for testname,testclip in testpolys_degens.items():
        for mode in ("intersect","union","difference"):
            test_draw(testname, subjpoly, testclip, mode)
    for testname,testclip in testpolys_nextto_almostsame.items():
        for mode in ("intersect","union","difference"):
            test_draw(testname, subjpoly, testclip, mode)
