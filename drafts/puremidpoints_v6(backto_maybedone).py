# -*- coding: UTF-8 -*-
# Efficient Clipping of Arbitrary Polygons
#
# Copyright (c) 2011, 2012 Helder Correia <helder.mc@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""
# Greiner-Hormann Polygon Clipping with support for degenerates

This is a fork aimed to improve Helder Correia's pure-Python Greiner-Hormann implementation for polygon clipping. Partly for educational purposes and partly for portable pure-Python clipping.

Status: Incomplete/unstable.

Fork author: Karim Bahgat <karim.bahgat.norway@gmail.com>

-----------------------------------------------------------

# Efficient Clipping of Arbitrary Polygons

Based on the paper "Efficient Clipping of Arbitrary Polygons" by Günther
Greiner (greiner[at]informatik.uni-erlangen.de) and Kai Hormann
(hormann[at]informatik.tu-clausthal.de), ACM Transactions on Graphics
1998;17(2):71-83.

Available at: http://www.inf.usi.ch/hormann/papers/Greiner.1998.ECO.pdf

You should have received the README file along with this program.
If not, see <https://github.com/helderco/polyclip>
"""


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
        self.alpha = alpha          # intersection point's relative distance from previous vertex
        self.intersect = intersect  # True if vertex is an intersection
        self.degen = degen
        self.checked = checked      # True if the vertex has been checked (last phase)

    def isInside(self, poly):
        if testLocation(self, poly) in ("in","on"):
            return True
        else: return False

    def setChecked(self):
        self.checked = True
        if self.neighbour and not self.neighbour.checked:
            self.neighbour.setChecked()

    def __repr__(self):
        """String representation of the vertex for debugging purposes."""
        return "(%.2f, %.2f) <-> %s(%.2f, %.2f)%s <-> (%.2f, %.2f) %s" % (
            self.prev.x, self.prev.y,
            'i' if self.intersect else ' ',
            self.x, self.y,
            ('e' if self.entry else 'x') if self.intersect else ' ',
            self.next.x, self.next.y,
            ' !' if self.intersect and not self.checked else ''
            )


class Polygon(object):
    """Manages a circular doubly linked list of Vertex objects that represents a polygon."""

    first = None

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

    def next(self, v):
        """Return the next non intersecting vertex after the one specified."""
        c = v
        while c.intersect:
            c = c.next
        return c

    @property
    def first_intersect(self):
        """Return the first unchecked intersection point in the polygon."""
        for v in self.iter():
            if v.intersect and not v.checked:
                break
        return v

    @property
    def points(self):
        """Return the polygon's points as a list of tuples (ordered coordinates pair)."""
        p = []
        for v in self.iter():
            p.append((v.x, v.y))
        return p

    def unprocessed(self):
        """Check if any unchecked intersections remain in the polygon."""
        for v in self.iter():
            if v.intersect and not v.checked:
                yield True

    def union(self, clip):
        return self.clip(clip, False, False)

    def intersect(self, clip):
        return self.clip(clip, True, True)

    def difference(self, clip):
        return self.clip(clip, False, True)

    def clip(self, clip, s_entry, c_entry):
        """Clip this polygon using another one as a clipper.

        This is where the algorithm is executed. It allows you to make
        a UNION, INTERSECT or DIFFERENCE operation between two polygons.

        Given two polygons A, B the following operations may be performed:

        A|B ... A OR B  (Union of A and B)
        A&B ... A AND B (Intersection of A and B)
        A\B ... A - B
        B\A ... B - A

        The entry records store the direction the algorithm should take when
        it arrives at that entry point in an intersection. Depending on the
        operation requested, the direction is set as follows for entry points
        (f=forward, b=backward; exit points are always set to the opposite):

              Entry
              A   B
              -----
        A|B   b   b
        A&B   f   f
        A\B   b   f
        B\A   f   b

        f = True, b = False when stored in the entry record
        """

        #   CHANGES TO ORIGINAL ALGORITHM:
        #       during intersection phase mark .on flag as True if alphas are 0 or 1
        #       then mark entry/exit as usual, but if .on flag then check prev and next and neighbour location first


        # detect clip mode
        unionmode = not s_entry and not c_entry
        intersectionmode = s_entry and c_entry
        differencemode = not s_entry and c_entry

        # prep by removing repeat of startpoint at end
        first = self.first
        last = first.prev
        if last.x == first.x and last.y == first.y:
            first.prev = last.prev
            last.prev.next = first
        first = clip.first
        last = first.prev
        if last.x == first.x and last.y == first.y:
            first.prev = last.prev
            last.prev.next = first
        
        # phase one - find intersections
        # ------------------------------
        anyintersection = False
        for s in self.iter(): # for each vertex Si of subject polygon do
            if not s.intersect:
                for c in clip.iter(): # for each vertex Cj of clip polygon do
                    if not c.intersect:
                        try:
                            
                            i, alphaS, alphaC = intersect_or_on(s, self.next(s.next),
                                                          c, clip.next(c.next))

                            s_between = (0 < alphaS < 1)
                            c_between = (0 < alphaC < 1)
                            if s_between and c_between:
                                #both subj and clip intersect each other somewhere in the middle
                                iS = Vertex(i, alphaS, intersect=True, entry=False)
                                iC = Vertex(i, alphaC, intersect=True, entry=False)
                                self.insert(iS, s, self.next(s.next))
                                clip.insert(iC, c, clip.next(c.next))
                            else:
                                if s_between:
                                    #subj line is touched by the start or stop point of a line from the clip polygon, so insert and mark that intersection as a degenerate
                                    iS = Vertex(i, alphaS, intersect=True, entry=False, degen=True)
                                    self.insert(iS, s, self.next(s.next))
                                elif alphaS == 0:
                                    #subj line starts at intersection, so mark the "degen"-flag, and replace vertex instead of inserting
                                    iS = Vertex(i, alphaS, intersect=True, entry=False, degen=True)
                                    self.replace(s, iS)
                                elif alphaS == 1:
                                    #subj line ends at intersection, so mark the "degen"-flag, and replace vertex instead of inserting
                                    iS = Vertex(i, alphaS, intersect=True, entry=False, degen=True)
                                    self.replace(self.next(s.next), iS)
                                if c_between:
                                    #clip line is touched by the start or stop point of a line from the subj polygon, so insert and mark that intersection as a degenerate
                                    iC = Vertex(i, alphaC, intersect=True, entry=False, degen=True)
                                    clip.insert(iC, c, clip.next(c.next))
                                elif alphaC == 0:
                                    #clip line starts at intersection, so mark the "degen"-flag, and replace vertex instead of inserting
                                    iC = Vertex(i, alphaC, intersect=True, entry=False, degen=True)
                                    clip.replace(c, iC)
                                elif alphaC == 1:
                                    #clip line ends at intersection, so mark the "degen"-flag, and replace vertex instead of inserting
                                    iC = Vertex(i, alphaC, intersect=True, entry=False, degen=True)
                                    clip.replace(clip.next(c.next), iC)
                                
                            iS.neighbour = iC
                            iC.neighbour = iS
                            
                            anyintersection = True
                            
                        except TypeError:
                            pass # this simply means intersect() returned None

        # phase one and a half - no intersections between subject and clip, so correctly return results
        # --------------------
        if not anyintersection: 
            resultpolys = []
            if unionmode: # union
                if clip.first.isInside(self):
                    # clip polygon is entirely inside subject, so just return subject shell
                    clipped = Polygon()
                    for s in self.iter():
                        clipped.add(Vertex(s))
                    polytuple = (clipped, [])
                    resultpolys.append(polytuple)
                elif self.first.isInside(clip):
                    # subject polygon is entirely inside clip, so just return clip shell
                    clipped = Polygon()
                    for c in clip.iter():
                        clipped.add(Vertex(c))
                    polytuple = (clipped, [])
                    resultpolys.append(polytuple)
                else:
                    #clip polygon is entirely outside subject, so return both
                    clipped = Polygon()
                    for s in self.iter():
                        clipped.add(Vertex(s))
                    polytuple = (clipped, [])
                    resultpolys.append(polytuple)
                    clipped = Polygon()
                    for c in clip.iter():
                        clipped.add(Vertex(c))
                    polytuple = (clipped, [])
                    resultpolys.append(polytuple)
            elif intersectionmode: # intersection
                if clip.first.isInside(self):
                    # clip polygon is entirely inside subject, so the intersection is only the clip polygon
                    clipped = Polygon()
                    for c in clip.iter():
                        clipped.add(Vertex(c))
                    polytuple = (clipped, [])
                    resultpolys.append(polytuple)
                elif self.first.isInside(clip):
                    # subject polygon is entirely inside clip, so the intersection is only the subject polygon
                    clipped = Polygon()
                    for s in self.iter():
                        clipped.add(Vertex(s))
                    polytuple = (clipped, [])
                    resultpolys.append(polytuple)
                else:
                    #clip polygon is entirely outside subject, so no intersection to return
                    pass
            elif differencemode: # difference
                if clip.first.isInside(self):
                    # clip polygon is entirely inside subject, so the difference is subject with clip as a hole
                    clipped = Polygon()
                    for s in self.iter():
                        clipped.add(Vertex(s))
                    hole = Polygon()
                    for c in clip.iter():
                        hole.add(Vertex(c))
                    polytuple = (clipped, [hole])
                    resultpolys.append(polytuple)
                elif self.first.isInside(clip):
                    # subject polygon is entirely inside clip, so there is no difference
                    pass
                else:
                    #clip polygon is entirely outside subject, so difference is simply the subject
                    clipped = Polygon()
                    for s in self.iter():
                        clipped.add(Vertex(s))
                    polytuple = (clipped, [])
                    resultpolys.append(polytuple)
            # no need to continue so just return result
            return resultpolys







        # phase two - identify entry/exit points
        # --------------------------------------

        #firstloc = testLocation(self.first, clip)
        #s_entry ^= firstloc in ("in","on")
        print "starts as ",s_entry
        for s in self.iter():
            if s.intersect:
                print "intersection"
                print "\t",s
                # intersection is degenerate, is the start/endpoint of a line
                # so maybe delete intersection flag based on prev/next locations
                prevloc = testLocation(s.prev, clip)
                nextloc = testLocation(s.next, clip)
                if prevloc == "on" or nextloc == "on":
                    prevmid = Vertex(((s.x+s.prev.x)/2.0,(s.y+s.prev.y)/2.0))
                    prevloc = testLocation(prevmid, clip)
                    nextmid = Vertex(((s.x+s.next.x)/2.0,(s.y+s.next.y)/2.0))
                    nextloc = testLocation(nextmid, clip)
                print "\t %s -> degenintsec -> %s" %(prevloc,nextloc)
                if prevloc == "out":
                    if nextloc == "out":
                        #just touching
                        s.intersect = False
                    elif nextloc == "in":
                        s.entry = s_entry
                    elif nextloc == "on":
                        s.entry = s_entry
                elif prevloc == "in":
                    #union and difference should never go inside the other polygon
                    #so this should only happen for intersectmode...
                    if nextloc == "in":
                        #just touching
                        s.intersect = False
                    elif nextloc == "out":
                        s.entry = not s_entry
                    elif nextloc == "on":
                        s.entry = not s_entry
                elif prevloc == "on":
                    if nextloc == "on":
                        s.intersect = False
                    elif nextloc == "out":
                        prev = s.prev
                        while not prev.intersect:
                            prev = prev.prev
                        if prev.entry != s_entry:
                            s.intersect = False
                        else:
                            s.entry = not s_entry
                    elif nextloc == "in":
                        prev = s.prev
                        while not prev.intersect:
                            prev = prev.prev
                        if prev.entry == s_entry:
                            s.intersect = False
                        else:
                            s.entry = s_entry
                    #elif nextloc == "out":# and s_entry:
                    #    s_entry = False
##                            #looking to enter (ie came from outside) but went back out again
##                            if s_entry: 
##                                if unionmode:
##                                    s.intersect = False
##                            #looking to exit (ie came from inside) and indeed went out
##                            elif not s_entry:
##                                if unionmode or differencemode:
##                                    s.intersect = False
                    #elif nextloc == "in": #and not s_entry:
                    #    s_entry = True
##                            #looking to enter (ie came from outside) and indeed went in
##                            if s_entry: 
##                                if unionmode:
##                                    s.intersect = False
##                            #looking to exit (ie came from inside) but went back inside again
##                            elif not s_entry:
##                                s.intersect = False
                #if intersect wasnt deleted, set and toggle entry flag
                if s.intersect:
                    print "\t entering = ", s.entry
            else:
                print "vertex"
                print "\t",s
        # then do same for clip polygon
        print "---"
        #firstloc = testLocation(clip.first, self)
        #c_entry ^= firstloc in ("in","on")
        print "starts as ",c_entry
        for c in clip.iter():
            if c.intersect:
                print "intersection"
                print "\t",c
                # intersection is degenerate, is the start/endpoint of a line
                # so maybe delete intersection flag based on prev/next locations
                prevloc = testLocation(c.prev, self)
                nextloc = testLocation(c.next, self)
                if prevloc == "on" or nextloc == "on":
                    prevmid = Vertex(((c.x+c.prev.x)/2.0,(c.y+c.prev.y)/2.0))
                    prevloc = testLocation(prevmid, self)
                    nextmid = Vertex(((c.x+c.next.x)/2.0,(c.y+c.next.y)/2.0))
                    nextloc = testLocation(nextmid, self)
                print "\t %s -> degenintsec -> %s" %(prevloc,nextloc)
                if prevloc == "out":
                    if nextloc == "out":
                        #just touching
                        c.intersect = False
                    elif nextloc == "in":
                        c.entry = c_entry
                    elif nextloc == "on":
                        c.entry = c_entry
                elif prevloc == "in":
                    #union and difference should never go inside the other polygon
                    #so this should only happen for intersectmode...
                    if nextloc == "in":
                        #just touching
                        c.intersect = False
                    elif nextloc == "out":
                        c.entry = not c_entry
                    elif nextloc == "on":
                        c.entry = not c_entry
                elif prevloc == "on":
                    if nextloc == "on":
                        c.intersect = False
                    elif nextloc == "out":
                        prev = c.prev
                        while not prev.intersect:
                            prev = prev.prev
                        if prev.entry != c_entry:
                            c.intersect = False
                        else:
                            c.entry = not c_entry
                    elif nextloc == "in":
                        prev = c.prev
                        while not prev.intersect:
                            prev = prev.prev
                        if prev.entry == c_entry:
                            c.intersect = False
                        else:
                            c.entry = c_entry
                    #elif nextloc == "out":# and s_entry:
                    #    s_entry = False
##                            #looking to enter (ie came from outside) but went back out again
##                            if s_entry: 
##                                if unionmode:
##                                    s.intersect = False
##                            #looking to exit (ie came from inside) and indeed went out
##                            elif not s_entry:
##                                if unionmode or differencemode:
##                                    s.intersect = False
                    #elif nextloc == "in": #and not s_entry:
                    #    s_entry = True
##                            #looking to enter (ie came from outside) and indeed went in
##                            if s_entry: 
##                                if unionmode:
##                                    s.intersect = False
##                            #looking to exit (ie came from inside) but went back inside again
##                            elif not s_entry:
##                                s.intersect = False
                #if intersect wasnt deleted, set and toggle entry flag
                if c.intersect:
                    print "\t entering = ", c.entry
                    # check if neighbours have same entry flag
                    #if c.entry == c.neighbour.entry:
                    #    c.neighbour.intersect = False
            else:
                print "vertex"
                print "\t",c







        # phase three - construct a list of clipped polygons
        # --------------------------------------------------
        resultpolys = []
        for _ in self.unprocessed():
            print "traverse new outpolygon"
            current = self.first_intersect
            clipped = Polygon()
            clipped.add(Vertex(current))
            while True:
                current.setChecked()
                print "flag is %s for %s" %(current.entry, current)
                if current.entry:
                    while True:
                        current = current.next
                        print "next",current
                        clipped.add(Vertex(current))
                        if current.intersect:
                            print "reached intersection"
                            break
                else:
                    while True:
                        current = current.prev
                        print "prev",current
                        clipped.add(Vertex(current))
                        if current.intersect:
                            print "reached intersection"
                            break

                current = current.neighbour
                if current.checked:
                    break

            polytuple = (clipped, [])
            resultpolys.append(polytuple)

        # finally, sort into exteriors and holes
        for pindex,(polyext,polyholes) in enumerate(resultpolys):
            for otherext,otherholes in resultpolys:
                if polyext == otherext:
                    continue # don't compare to self
                if polyext.first.isInside(otherext):
                    otherholes.append(polyext) #poly is within other so make into a hole
                    del resultpolys[pindex] #and delete poly from being an independent poly
        return resultpolys

    def __repr__(self):
        """String representation of the polygon for debugging purposes."""
        count, out = 1, "\n"
        for s in self.iter():
            out += "%02d: %s\n" % (count, str(s))
            count += 1
        return out

    def iter(self):
        """Iterator generator for this doubly linked list."""
        s = self.first
        while True:
            yield s
            s = s.next
            if s == self.first:
                return




def intersect_or_on(s1, s2, c1, c2):
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
    w =0
    for v in polygon.iter():
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

def clip_polygon(subject, clipper, operation = 'difference'):
    """
    Higher level function for clipping two polygons (from a list of points).
    Since input polygons are lists of points, output is also in list format.
    Each polygon in the resultlist is a tuple of: (polygon exterior, list of polygon holes)
    """
    Subject = Polygon()
    Clipper = Polygon()

    for s in subject:
        Subject.add(Vertex(s))

    for c in clipper:
        Clipper.add(Vertex(c))

    clipped = Clipper.difference(Subject)\
    if operation == 'reversed-diff'\
    else Subject.__getattribute__(operation)(Clipper)

    clipped = [(ext.points,[hole.points for hole in holes]) for ext,holes in clipped]
    return clipped





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
    
    def test_draw(testname, subjpoly, clippoly, mode):
        t = time.time()
        resultpolys = clip_polygon(subjpoly,clippoly,mode)
        print "finished:",resultpolys,time.time()-t
        crs = pydraw.CoordinateSystem([-1,-1,11,11])
        img = pydraw.Image(400,400, crs=crs)
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
    
