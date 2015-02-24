# Martinez et al polygon clipping
# Much faster than both Vatti and Greiner for large samples,
# ...and only barely slower for very small samples
# http://www.cs.ucr.edu/~vbz/cs230papers/martinez_boolean.pdf
# example implementation: https://github.com/akavel/polyclip-go


# Classes and helpers

def pairwise(points):
    a = (p for p in points)
    b = (p for p in points)
    next(b)
    for curpoint,nextpoint in zip(a, b):
        yield curpoint, nextpoint


class Sweepline:
    def __init__(self):
        self.data = []

    def insert(self, edge):
        pass

    def next(self, edge):
        pass

    def prev(self, edge):
        pass
    

class Endpoint: # aka sweepevent
    def __init__(self, xy, other, polytype):
        self.x,self.y = xy
        self.other = other
        self.left = xy[0] < other[0]
        self.polytype = polytype
        self.inout = ""
        self.inside = ""
        self.edgetype = ""


def set_inside_flag(endpoint1, endpoint2):
    if not endpoint2:
        endpoint1.inside = endpoint1.inout = False
    elif endpoint1.polytype == endpoint2.polytype:
        endpoint1.inside = endpoint2.inside
        endpoint1.inout = endpoint2.inout
    else:
        endpoint1.inside = not endpoint2.inout
        endpoint1.inout = endpoint2.inside


def possible_intersect(edge1, edge2):
    """Looks for intersections so that it can subdivide them
    and update the Q and S lists.
    """

    if edge1.polytype == edge2.polytype:
        # belong to same polytype, stop processing
        return
    
    intersect = intersect_test(edge1, edge2)
    
    if intersect:
        point,iS,iC = intersect
        if 0.0 < iS < 1.0 or 0.0 < iC < 1.0:
            # subdivide so the lines do not intersect...?
            # ie replace the isect point to the end of each edge
            # update Q and S
            pass

        else:
            # only intersect at one of their endpoints, stop processing
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




# The algorithm

def clip_polygons(subject, clip, type):

    # insert all edge endpoints into the priority queue    
    Q = []
    for start,end in pairwise(subject):
        Q.append( Endpoint(start, end, "subject") )
        Q.append( Endpoint(end, start, "subject") )
    for start,end in pairwise(clip):
        Q.append( Endpoint(start, end, "clip") )
        Q.append( Endpoint(end, start, "clip") )
        
    # sort by x asc, then y desc
    Q = sorted(Q, key=lambda EP: EP.y, reverse=True)
    Q = sorted(Q, key=lambda EP: EP.x)

    # create the sweepline
    S = Sweepline()

    # loop
    raw = []
    while Q:
        endpoint = Q.pop(0)
        if endpoint.left: # left endpoint
            # some lists
            edge = S.insert(endpoint)
            set_inside_flag(S.prev(edge))
            # intersections
            possible_intersect(edge, S.next (edge) )
            possible_intersect(edge, S.prev (edge) )

        else: # right endpoint
            # move forward
            edge = S.find(endpoint.other)
            next = S.next(edge)
            prev = S.prev(edge)
            # some adding
            if endpoint.inside:
                raw.extend(endpoint.segment) # intersection
            elif not endpoint.inside:
                raw.extend(endpoint.segment) # union
            # some cleanup
            S.erase(edge)
            possible_intersect(prev, next)

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
            chain = edge
            chains.append(chain)

                






    

