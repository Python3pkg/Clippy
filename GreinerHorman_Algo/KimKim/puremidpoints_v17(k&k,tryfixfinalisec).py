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

# FINAL BEST IDEA, IMPLEMENTED NOW BUT DOESNT WORK CUS INTERSECTION STAGE NOT FINDING ALL INTERSECTIONS
# USE PREV AND NEXT MIDPOINT LOCS FOR DETERMINING ENTRY FLAG
# NORMAL RULES, EXCEPT FOR INTERSECTIONMODE TURN OFF INTERSECTIONFLAGS FOR OUT-ON-ON and ON-ON-OUT BC THEY ARE JUST TANGENT AND NOT RELATED TO INSIDES
# FINALLY WHEN TRAVERSING, AFTER COMPLETING ONE POLY, SEARCH FOR NEXT ISECT THAT IS UNCHECK IN BOTH CURRENT AND NEIGHRBOUR

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



DEBUG = False


class Vertex(object):
    """Node in a circular doubly linked list.

    This class is almost exactly as described in the paper by Günther/Greiner.
    """

    def __init__(self, vertex, alpha=0.0, intersect=False, entry=None, checked=False, degen=False, orig=True):
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
        self.checked = checked      # True if the vertex has been checked (last phase)
        self.couple = None
        self.cross_change = None
        self.orig = orig

    @property
    def xy(self):
        return self.x, self.y

    def isInside(self, poly):
        if testLocation(self, poly) in ("in","on"):
            return True
        else: return False

    def setChecked(self):
        self.checked = True
        if self.neighbour and not self.neighbour.checked:
            self.neighbour.setChecked()

    def copy(self):
        copy = Vertex(self)     # point coordinates of the vertex
        copy.next = self.__next__        # reference to the next vertex of the polygon
        copy.prev = self.prev            # reference to the previous vertex of the polygon
        copy.neighbour = self.neighbour       # reference to the corresponding intersection vertex in the other polygon
        copy.entry = self.entry          # True if intersection is an entry point, False if exit
        copy.alpha = self.alpha          # intersection point's relative distance from previous vertex
        copy.intersect = self.intersect  # True if vertex is an intersection
        copy.couple = self.couple
        copy.cross_change = self.cross_change
        copy.checked = self.checked
        return copy
    
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
        # when replacing old normal vertice with new intersection vertice at same xy
        # only changes the attributes in place
        old.intersect = new.intersect
        old.x,old.y = new.x,new.y
        old.neighbour = new.neighbour
        old.neighbour.neighbour = old
        old.entry = new.entry
        old.alpha = new.alpha
        
##        new.next = old.next
##        new.prev = old.prev
##        if old == self.first:
##            #print "replaced first", self.first, new
##            self.first = new
##        old.prev.next = new
##        old.next.prev = new

    def insert(self, vertex, start, end):
        """Insert and sort a vertex between a specified pair of vertices.

        This function inserts a vertex (most likely an intersection point)
        between two other vertices (start and end). These other vertices
        cannot be intersections (that is, they must be actual vertices of
        the original polygon). If there are multiple intersection points
        between the two vertices, then the new vertex is inserted based on
        its alpha value.
        """
        if vertex.xy == start.xy:
            copy = vertex.copy()
            self.replace(start, copy)
            return # dont process further

        elif vertex.xy == end.xy:
            copy = vertex.copy()
            self.replace(end, copy)
            return # dont process further

        # position based on alpha
        curr = start
        while curr != end and curr.alpha < vertex.alpha:
            curr = curr.__next__

        if vertex.xy == curr.prev.xy:
##            if vertex.xy == curr.xy: self.replace(curr, vertex)
##            elif vertex.xy == curr.prev.xy: self.replace(curr, vertex.prev)
            vertex.neighbour.neighbour = curr.prev
            return # dont do it if same as a previously inserted intersection

        if vertex.xy == curr.xy:
##            if vertex.xy == curr.xy: self.replace(curr, vertex)
##            elif vertex.xy == curr.prev.xy: self.replace(curr, vertex.prev)
            vertex.neighbour.neighbour = curr
            return # dont do it if same as a previously inserted intersection
        
        vertex.next = curr
        vertex.prev = curr.prev
        vertex.next.prev = vertex
        vertex.prev.next = vertex
        #print "inserted",vertex

    def next_orig(self, v):
        """Return the next original vertex after the one specified."""
        c = v.__next__
        while not c.orig:
            c = c.__next__
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

        # TODO: maybe also remove repeat points anywhere?
        # ...
        
        # phase one - find intersections
        # ------------------------------
        anyintersection = False
        
        s = self.first
        while s != self.first.prev:
            c = clip.first
            while c != clip.first.prev:
                s_next_orig = self.next_orig(s)
                c_next_orig = clip.next_orig(c)
                res = intersect_or_on(s, s_next_orig,
                                        c, c_next_orig)
                if res:
                    i, alphaS, alphaC = res
                    
##                    if alphaS == 0:
##                        iS = s
##                        iS.intersect = True
##                        iS.entry = False
##                    elif alphaS == 1:
##                        iS = s_next_orig
##                        iS.intersect = True
##                        iS.entry = False
##                    else:
##                        # insert
##                        iS = Vertex(i, alphaS, intersect=True, entry=False, orig=False)
##                        v = s.next
##                        while v != s_next_orig and v.alpha > alphaS:
##                            v = v.next
##                        print "insert S",s.xy,iS.xy,v.xy
##                        s.next = iS
##                        iS.prev = s
##                        v.prev = iS
##                        iS.next = v
##
##                    if alphaC == 0:
##                        iC = c
##                        iC.intersect = True
##                        iC.entry = False
##                    elif alphaC == 1:
##                        iC = c_next_orig
##                        iC.intersect = True
##                        iC.entry = False
##                    else:
##                        # insert
##                        iC = Vertex(i, alphaC, intersect=True, entry=False, orig=False)
##                        v = c.next
##                        while v != c_next_orig and v.alpha > alphaC:
##                            v = v.next
##                        print "insert C",c.xy,iC.xy,v.xy
##                        c.next = iC
##                        iC.prev = c
##                        v.prev = iC
##                        iC.next = v

                    iS = Vertex(i, alphaS, intersect=True, entry=False, orig=False)
                    iC = Vertex(i, alphaC, intersect=True, entry=False, orig=False)

                    iS.neighbour = iC
                    iC.neighbour = iS

                    self.insert(iS, s, s_next_orig)
                    print("insert S",s.xy,iS.xy,s_next_orig.xy)
                    clip.insert(iC, c, c_next_orig)
                    print("insert C",c.xy,iC.xy,c_next_orig.xy)
                    
                    anyintersection = True
                
                c = clip.next_orig(c)
            s = self.next_orig(s)

##        s_intsecs = []
##        c_intsecs = []
##        for s in self.iter(): # for each vertex Si of subject polygon do
##            for c in clip.iter(): # for each vertex Cj of clip polygon do
##                try:
##                    #print "find isect %s - %s and %s - %s" %(s.xy, self.next(s.next).xy, c.xy, clip.next(c.next).xy )
##                    s_next = self.next_nonintsec(s)
##                    c_next = clip.next_nonintsec(c)
##                    i, alphaS, alphaC = intersect_or_on(s, s_next,
##                                                        c, c_next)
##                    
##                    iS = Vertex(i, alphaS, intersect=True, entry=False)
##                    iC = Vertex(i, alphaC, intersect=True, entry=False)
##
##                    iS.neighbour = iC
##                    iC.neighbour = iS
##                    
##                    s_intsecs.append( (iS, alphaS, s, s_next) )
##                    c_intsecs.append( (iC, alphaC, c, c_next) )
##     
##                    anyintersection = True
##                    
##                except TypeError:
##                    pass # this simply means intersect() returned None
##
##        # insert intersections into originals
##        for iS,a,s,s_next in reversed(s_intsecs):
##            if a == 0:
##                self.replace(s, iS)
##            elif a == 1:
##                self.replace(s_next, iS)
##            else:
##                self.insert(iS, s, s_next)
##        for iC,a,c,c_next in reversed(c_intsecs):
##            if a == 0:
##                self.replace(c, iC)
##            elif a == 1:
##                self.replace(c_next, iC)
##            else:
##                clip.insert(iC, c, c_next)

        #print "testing if insert was done correctly"
        for s in self.iter():
            #print s
            pass
        #print "and"
        for c in clip.iter():
            #print c
            pass
                    

        # phase one and a half - no intersections between subject and clip, so correctly return results
        # --------------------
        def specialcase_insidetest():
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
        
        if not anyintersection: 
            return specialcase_insidetest()







        # phase two - identify entry/exit points
        # --------------------------------------

        # From K&K
        
        def mark_flags(poly, c, c_entry):
            "c and c_entry are not actually the clip, can be for both s and c, just too lazy to change."
            #print "intersection"
            #print "\t",c
            # intersection is degenerate, is the start/endpoint of a line
            # so maybe delete intersection flag based on prev/next locations
            prevloc = testLocation(c.prev, poly)
            nextloc = testLocation(c.__next__, poly)
            if prevloc == "on" or nextloc == "on":
                prevmid = Vertex(((c.x+c.prev.x)/2.0,(c.y+c.prev.y)/2.0))
                prevloc = testLocation(prevmid, poly)
                nextmid = Vertex(((c.x+c.next.x)/2.0,(c.y+c.next.y)/2.0))
                nextloc = testLocation(nextmid, poly)
            if prevloc == "in" or nextloc == "in":
                poly.anyinside = True
            #print "\t %s -> degenintsec -> %s" %(prevloc,nextloc)
            if prevloc == "out":
                if nextloc == "out":
                    #just touching
                    c.entry = "en/ex" if c_entry else "ex/en"
                elif nextloc == "in":
                    c.entry = "en" if c_entry else "ex"
                elif nextloc == "on":
                    c.entry = "en" if c_entry else "ex"
            elif prevloc == "in":
                #union and difference should never go inside the other polygon
                #so this should only happen for intersectmode...
                if nextloc == "in":
                    #just touching
                    c.entry = "ex/en" if c_entry else "en/ex"
                elif nextloc == "out":
                    c.entry = "ex" if c_entry else "en"
                elif nextloc == "on":
                    c.entry = "ex" if c_entry else "en"
            elif prevloc == "on":
                if nextloc == "on":
                    c.entry = None
                elif nextloc == "out":
                    c.entry = "ex" if c_entry else "en"
                elif nextloc == "in":
                    c.entry = "en" if c_entry else "ex"

        self.anyinside = False

        # set clip
        prevsingle = None
        for c in clip.iter():
            if c.intersect:
                mark_flags(self, c, c_entry)
                # set couple
                if c.entry in ("ex","en"):
                    if prevsingle and c.entry == prevsingle.entry:
                        c.couple = prevsingle
                        prevsingle.couple = c
                    prevsingle = c
                # set crosschange
                # some modifications based on implementation in Qt clipper source code
                #if c.entry == "en/ex" == c.neighbour.entry or c.entry == "ex/en" == c.neighbour.entry:
                if c.entry == "en/ex" or c.entry == "ex/en":
                    print("Maybe crosschange...")
                    # tri1
                    #a,b,c = c.neighbour.prev, c.prev, c.neighbour.next
                    a,b,c = c.neighbour.__next__, c.prev, c.neighbour.prev
                    dir1 = 0.5 * (a.x * (b.y-c.y) +
                                  b.x * (c.y-a.y) +
                                  c.x * (a.y-b.y))
                    # tri2
                    #a,b,c = c.neighbour.prev, c.prev, c.next
                    a,b,c = c.__next__, c.prev, c.neighbour.prev
                    dir2 = 0.5 * (a.x * (b.y-c.y) +
                                  b.x * (c.y-a.y) +
                                  c.x * (a.y-b.y))
                    print(dir1,dir2)
                    #if dir1 < 0 != dir2 < 0: # different orientation
                    if (dir1 * dir2) < 0: # different orientation means at least one negative, making the results less than 0
                        print("CROSSCHANGE!!!")
                        c.cross_change = True
                        c.neighbour.cross_change = True # not sure if should set neighbour too

        # maybe early abort
        if not self.anyinside and intersectionmode:
            return []

        # what about perfect overlap???
        # ...

        if False: #DEBUG:
            print("view clip entries")
            for c in clip.iter():
                print(c, c.entry)

        # find first isect where both neighbours have valid flag
        for c in clip.iter():
            if c.entry:
                s = c.neighbour
                mark_flags(clip, s, s_entry)
                if s.entry:
                    first_c = c
                    first_s = s
                    # print 777,s.entry
                    break

        else:
            return specialcase_insidetest()
            #raise Exception("weird special case, no neighbours that both have flag left")
        
        # autoset subj, if neighbour of first is different, then set all as opposite
        # TODO: how deal with s_entry in case of different modes...?
        print("view first")
        print(first_c, first_c.entry)
        print(first_s, first_s.entry)
        if first_c.entry != first_s.entry: # and s_entry: # this is the behaviour for standard intersect mode, otherwise flip, hence the s_entry
            for c in clip.iter():
                if c.entry:
                    if c.entry == "en": c.neighbour.entry = "ex"
                    elif c.entry == "ex": c.neighbour.entry = "en"
                    elif c.entry == "en/ex": c.neighbour.entry = "ex/en"
                    elif c.entry == "ex/en": c.neighbour.entry = "en/ex"
                    
        # else set all same
        else:
            for c in clip.iter():
                if c.entry:
                    c.neighbour.entry = c.entry

        # set couple for subj (not sure if needed)
        prevsingle = None
        for s in self.iter():
            if s.entry:
                if s.entry in ("ex","en"):
                    if prevsingle and s.entry == prevsingle.entry:
                        s.couple = prevsingle
                        prevsingle.couple = s
                    prevsingle = s

        if False: #DEBUG:       
            print("view subj entries")
            for s in self.iter():
                print(s, s.entry)





        # phase three - construct a list of clipped polygons
        # --------------------------------------------------

        ######
        # Defs
        def next_unprocessed(vert):
            origvert = vert
            while vert:
                if vert.entry and not (vert.checked or vert.neighbour.checked):
                    #print "vert, found next unproc", vert, vert.checked, vert.neighbour.checked
                    if vert.couple:
                        # rule 1
                        if vert.couple.entry and vert.entry:
                            # rule 2
                            if vert.couple.entry == "en" and vert.entry == "en":
                                return vert.couple
                            elif vert.couple.entry == "ex" and vert.entry == "ex":
                                return vert
                    # rule 3
                    else:
                        return vert
                
                vert = vert.__next__
                
                if vert == origvert:
                    # if returned to first, return None
                    return None

        def DeleteFlag1(cur, stat):
            if cur.entry == "en/ex":
                cur.entry = None
                if cur.cross_change:
                    if stat == "D3":
                        return "D3"
                    else:
                        return "D4"
                if stat == "D3":
                    return "D4"
                else:
                    return "D3"
            if cur.entry == "ex/en":
                if stat == "D3":
                    cur.entry = "en"
                    return "D2"
                else:
                    cur.entry = "ex"
                    return "D1"
            if cur.entry == "en":
                cur.entry = None
                return "D1"
            if cur.entry == "ex":
                cur.entry = None
                return "D2"

        def DeleteFlag2(cur, prev, stat):
            if cur.entry == "en/ex":
                if stat == "D1":
                    cur.entry = "ex"
                else:
                    cur.entry = "en"
                if cur.cross_change:
                    if stat == "D1":
                        return "D4"
                    else:
                        return "D3"
                if stat == "D1":
                    return "D3"
                else:
                    return "D4"
            if cur.entry == "ex/en":
                if stat == "D1":
                    cur.entry = "en"
                else:
                    cur.entry = "ex"
                if cur.cross_change:
                    if stat == "D1":
                        return "D4"
                    else:
                        return "D3"
                if stat == "D1":
                    return "D3"
                else:
                    return "D4"
            if cur.entry == "en":
                cur.entry = None
                if stat == "D1" and cur.couple and prev.couple == cur:
                    return "D1"
                if stat == "D1":
                    return "D3"
                else:
                    return "D4"
            if cur.entry == "ex":
                cur.entry = None
                if stat != "D1" and cur.couple and prev.couple == cur:
                    return "D2"
                else:
                    if stat == "D1":
                        return "D3"
                    else:
                        return "D4"

        def proceed(cur, stat):
            cur.checked = True
            if stat == "D1":
                clipped.add(Vertex(cur))
                return cur.__next__
            elif stat == "D2":
                clipped.add(Vertex(cur))
                return cur.prev
            else:
                return cur.neighbour

        ####
        resultpolys = []

        self.first.checked = True
        cur = prev = start = next_unprocessed(self.first)

        while cur:
            # each new polygon
            print("new poly")

            stat = DeleteFlag1(cur, "D3")
            if DEBUG: print("v", cur, cur.entry, stat)
            clipped = Polygon()
            cur = proceed(cur, stat)

            # collect vertexes
            while cur != start:
                if DEBUG: print("v", cur, cur.entry, stat)
                if cur.entry:
                    if stat == "D1" or stat == "D2":
                        stat = DeleteFlag2(cur, prev, stat)
                    else:
                        stat = DeleteFlag1(cur, stat)
                    prev = cur
                cur = proceed(cur, stat)

            # return to first vertex
            clipped.add(Vertex(clipped.first))

            print(clipped)

            resultpolys.append((clipped,[]))
            cur = prev = start = next_unprocessed(self.first)



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
            s = s.__next__
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

##    for s in Subject.iter():
##        if ("%.2f"%s.x == "12.14" and "%.2f"%s.y == "63.05") \
##        or ("%.2f"%s.x == "12.15" and "%.2f"%s.y == "63.00"):
##            print "PRE:"
##            print s
##            print s.neighbour
##            print "---"

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
##                        "smaller, outside":
##                        [(7,7),(7,9),(9,9),(9,7),(7,7)],
##                        "smaller, inside":
##                        [(2,2),(2,4),(4,4),(4,2),(2,2)],
##                        "larger, covering all":
##                        [(-1,-1),(-1,7),(7,7),(7,-1),(-1,-1)],
##                        "larger, outside":
##                        [(-10,-10),(-10,-70),(-70,-70),(-70,-10),(-10,-10)]
                        }

    # degenerate intersections
    testpolys_degens = {"degenerate, starts on edge intersection and goes inside":
                        [(0,5),(6,4),(10,4),(10,10),(4,10),(0,5)],
##                        "degenerate, starts on edge intersection and goes outside":
##                        [(5,6),(5.2,5.5),(5,5.4),(4.8,5.5)],
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

    DEBUG = False

    # test geo
    
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
##    test_draw("norway-sweden", norw, swed, "difference")
##
##    breakonpurpose

    # test basics

    def test_draw(testname, subjpoly, clippoly, mode):
        t = time.time()
        #print testname, mode
        resultpolys = clip_polygon(subjpoly,clippoly,mode)
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
    for testname,testclip in list(testpolys_normal.items()):
        print(testname)
        for mode in ("intersect","union","difference"):
            print(mode)
            test_draw(testname, subjpoly, testclip, mode)
    for testname,testclip in list(testpolys_degens.items()):
        print(testname)
        for mode in ("intersect","union","difference"):
            print(mode)
            test_draw(testname, subjpoly, testclip, mode)
    for testname,testclip in list(testpolys_nextto_almostsame.items()):
        print(testname)
        for mode in ("intersect","union","difference"):
            print(mode)
            test_draw(testname, subjpoly, testclip, mode)
    
