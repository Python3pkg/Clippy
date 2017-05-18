# Martinez et al polygon clipping
# Much faster than both Vatti and Greiner for large samples,
# ...and only barely slower for very small samples
# http://www.cs.ucr.edu/~vbz/cs230papers/martinez_boolean.pdf
# example implementation: https://github.com/akavel/polyclip-go

# STATUS: Mostly implemented but doesnt yet work even for simple cases
# not yet implemented support for overlapping edges or selfintersections

# Classes and helpers

def pairwise(points):
    a = (p for p in points)
    b = (p for p in points)
    next(b)
    for curpoint,nextpoint in zip(a, b):
        yield curpoint, nextpoint


class Sweepline:
    def __init__(self):
        self.edges = []

    def insert(self, edge):
        if len(self.edges) > 0:
            pos = 0
            for v in self.edges:
                if v.y > edge.y or (not edge.vertical and v.y == edge.y): # nonvertical are placed first for same points, otherwise after
                    break
                pos += 1
            self.edges.insert(pos,edge)
        else:
            self.edges.append(edge)

    def next(self, edge):
        for cur in self.edges:
            if cur.y > edge.y:
                break
        else:
            return None
        return cur

    def prev(self, edge):
        for cur in reversed(self.edges):
            if cur.y < edge.y:
                break
        else:
            return None
        return cur

    def erase(self, edge):
        self.edges.remove(edge)
    

class Endpoint: # aka sweepevent
    def __init__(self, xy, polytype):
        self.xy = xy
        self.x,self.y = xy
        self.other = None # must be set manually as a reference
        self.polytype = polytype
        self.inout = ""
        self.inside = ""
        self.edgetype = ""

    @property
    def vertical(self):
        return self.xy[0] == self.other.xy[0]

    @property
    def left(self):
        if self.vertical:
            return self.xy[1] < self.other.xy[1] # vertical edge, lower y is considered left
        else:
            return self.xy[0] < self.other.xy[0]

    def __str__(self):
        return ("Left" if self.left else "Right") + "Endpoint(%s, %s)" % self.xy + " --> (%s, %s)" % self.other.xy



def set_inside_flag(endpoint, prevendpoint):
    endpoint1, endpoint2 = endpoint, prevendpoint
    if not endpoint2:
        endpoint1.inside = endpoint1.inout = False
    elif endpoint1.polytype == endpoint2.polytype:
        endpoint1.inside = endpoint2.inside
        endpoint1.inout = endpoint2.inout
    else:
        endpoint1.inside = not endpoint2.inout
        endpoint1.inout = endpoint2.inside


def possible_subdivide(ep1, ep2, S, Q):
    """Looks for intersections so that it can subdivide them
    and update the Q and S lists.
    """

    if not (ep1 and ep2):
        # if only one edge, then no intersections
        return

    if ep1.polytype == ep2.polytype:
        # belong to same polytype, stop processing
        return
    
    intersect = intersect_or_on(ep1, ep2)
    
    if intersect:
        ipoint,alphaS,alphaC = intersect
        if 0.0 < alphaS < 1.0 or 0.0 < alphaC < 1.0:
            # subdivide so the lines do not intersect...?
            # ie replace the isect point to the end of each edge
            right = ep1.other
            ep1.other = Endpoint(ipoint, ep1.polytype)
            ep1.other.other = ep1
            assert ep1.left
            ep1ext = Endpoint(ipoint, ep1.polytype)
            ep1ext.other = right
            right.other = ep1ext
            assert ep1ext.left
            right = ep2.other
            ep2.other = Endpoint(ipoint, ep2.polytype)
            ep2.other.other = ep2
            assert ep2.left
            ep2ext = Endpoint(ipoint, ep2.polytype)
            ep2ext.other = right
            right.other = ep2ext
            assert ep2ext.left
            
            # update Q and S
            # TODO: The Q endpoints mut be inserted sorted bottom to top plus more
            
##            cur = Q[1] if len(Q) > 1 else Q[0]
##            curpos = 0
##            while cur.x == p2.x and cur.y < p2.y and not cur.left:
##                curpos += 1
##            Q.insert(pos-1,p2)
##
##            cur = Q[1] if len(Q) > 1 else Q[0]
##            curpos = 0
##            while cur.x == p4.x and cur.y < p4.y and not cur.left:
##                curpos += 1
##            Q.insert(pos-1,p4)
            
            for pos,v in enumerate(Q):
                if v.y > ep1ext.y or v.x != ep1ext.x or (v.xy == ep1ext.xy and not ep1ext.left and v.left):
                    assert pos >= 0
                    while v.xy == ep1ext.xy and ep1ext.left and v.left and S.edges.index(v) < S.edges.index(ep1ext):
                        pos += 1
                    Q.insert(pos,ep1ext)
                    break
            for pos,v in enumerate(Q):
                if v.y > ep2ext.y or v.x != ep2ext.x or (v.xy == ep2ext.xy and not ep2ext.left and v.left):
                    assert pos >= 0
                    while v.xy == ep2ext.xy and ep1ext.left and v.left and S.edges.index(v) < S.edges.index(ep2ext):
                        pos += 1
                    Q.insert(pos,ep2ext)
                    break
    
            S.insert(ep1ext)
            S.insert(ep2ext)

        else:
            # only intersect at one of their endpoints, stop processing
            return


def intersect_or_on(ep1, ep2):
    """Same as intersect(), except returns
    intersection even if degenerate.
    """
    s1, s2 = ep1, ep1.other
    c1, c2 = ep2, ep2.other
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




# The algorithm

def clip_polygons(subject, clip, type):

    # insert all edge endpoints into the priority queue    
    Q = []
    for start,end in pairwise(subject):
        ep1 = Endpoint(start, "subject")
        ep2 = Endpoint(end, "subject")
        ep1.other = ep2
        ep2.other = ep1
        Q.append( ep1 )
        Q.append( ep2 )
    for start,end in pairwise(clip):
        ep1 = Endpoint(start, "clip")
        ep2 = Endpoint(end, "clip")
        ep1.other = ep2
        ep2.other = ep1
        Q.append( ep1 )
        Q.append( ep2 )
        
    # sort by x asc, y asc (bottom to top), and right endpoints first
    Q = sorted(Q, key=lambda EP: (EP.x,EP.y,EP.left,EP.other.y))
    for q in Q:
        print("sorted",q)
        
    # create the sweepline
    S = Sweepline()

    # loop
    raw = []
    while Q:
        endpoint = Q.pop(0)
        print("processing",endpoint)
        if endpoint.left: # left endpoint
            # some lists
            S.insert(endpoint)
            set_inside_flag(endpoint, S.prev(endpoint))
            # intersections
            possible_subdivide(endpoint, S.next(endpoint), S, Q)
            possible_subdivide(endpoint, S.prev(endpoint), S, Q)

        else: # right endpoint
            # move forward
            #S.find(endpoint.other)
            next = S.next(endpoint.other)
            prev = S.prev(endpoint.other)
            #print prev,endpoint.other,next
            # some adding
            if endpoint.inside:
                raw.append(endpoint.xy) # intersection
                raw.append(endpoint.other.xy)
            elif not endpoint.inside:
                raw.append(endpoint.xy) # union
                raw.append(endpoint.other.xy)
            # some cleanup
            #print "remove",endpoint.other,"from",[str(e) for e in S.edges]
            #S.erase(endpoint) if endpoint in S.edges else S.erase(endpoint.other)
            if endpoint.other in S.edges: S.erase(endpoint.other) #right ends are processed before left, so could be the left hasnt been inserted yet
            possible_subdivide(prev, next, S, Q)

    # connect the final edges
    polys = []
    chains = []
    for edge in pairwise(raw):
        connect = []
        for chain in chains:
            if chain[0] in edge:
                connect.append(chain)
            elif chain[1] in edge:
                connect.append(chain)

        if len(connect) == 2:
            chain1,chain2 = connect
            chain1.extend(chain2)
            chain1.extend(edge)
            chains.remove(chain2)
            if chain1[0] == chain1[-1]:
                polys.append( chain1 )
                chains.remove(chain1)

        elif len(connect) == 1:
            chain = connect[0]
            chain.extend(edge)
            if chain[0] == chain[-1]:
                polys.append( chain )
                chains.remove(chain)

        else:
            chain = list(edge)
            chains.append(chain)

    return [(p,[]) for p in polys]




                


if __name__ == "__main__":
    """
    Test and visualize various polygon overlap scenarios.
    Visualization requires the pure-Python PyDraw library from
    https://github.com/karimbahgat/PyDraw
    """
    
    subjpoly = [(0,0),(6,0),(7,6),(1,6)]
    
    # normal intersections
    testpolys_normal = {"simple overlap":
                        [(0+4,0+4),(6+4,0+4),(7+4,6+4),(1+4,6+4)],
                        "jigzaw overlap":
                        [(1,4),(3,8),(5,4),(6,10),(2,10)],
##                        "smaller, outside":
##                        [(7,7),(7,9),(9,9),(9,7),(7,7)],
##                        "smaller, inside":
##                        [(2,2),(2,4),(4,4),(4,2),(2,2)],
##                        "larger, covering all":
##                        [(-1,-1),(-1,7),(7,7),(7,-1),(-1,-1)],
##                        "larger, outside":
##                        [(-10,-10),(-10,-70),(-70,-70),(-70,-10),(-10,-10)]
                        }

    #run operation

    import os
    import time
    import pydraw

    def test_draw(testname, subjpoly, clippoly, mode):
        t = time.time()
        print(testname, mode)
        resultpolys = clip_polygons(subjpoly,clippoly,mode)
        print("finished:",resultpolys,time.time()-t)
        crs = pydraw.CoordinateSystem([-1,-1,11,11])
        img = pydraw.Image(300,300, crs=crs)
        img.drawpolygon(subjpoly, fillcolor=(222,0,0))
        img.drawpolygon(clippoly, fillcolor=(0,222,0))
        for ext,holes in resultpolys:
            img.drawpolygon(ext,holes)
        img.drawgridticks(1,1)
        img.view() #save("test_output/"+testname+"-"+mode+".png")


    for testname,testclip in list(testpolys_normal.items()):
        print(testname)
        for mode in ("intersect","union","difference"):
            print(mode)
            test_draw(testname, subjpoly, testclip, mode)

    

