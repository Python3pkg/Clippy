# Greiner-Hormann Polygon Clipping with support for degenerates
Fork Author: Karim Bahgat

This is a fork aimed to improve Helder Correia's pure-Python Greiner-Hormann implementation for polygon clipping. Partly for educational purposes and partly for portable pure-Python clipping. 

## Improvement 1: Correctly handle totally inside and totally outside cases

When the subject and clip polygon are either totally inside or outside each other, the original implementation simply exited and returned an empty polygon (non-result) instead of returning the correct result based on the requested operation (eg a union of two totally separate polygons should return both).

__Status__: Complete

## Improvement 2: Differentiate exterior from holes in the result-polygons. 

Helder Correia's original implementation returned the result as a list of polygons, regardless if they were exteriors or holes, leaving it up to the user to differentiate them. In this fork each resultpolygon are divided into exterior and holes. 

__Status__: Almost complete

For now, holes are determined by testing the location of only the first vertex visavis the other polygon, but it could maybe be possible that the returned polygons cross each other, in which case this is not enough. Investigate if further hole-testing is required. 

## Improvement 3: Add support for degenerates

Since the original implementation did not support degenerates (when polygons share edges/points), the purpose of this fork is to add support for such degenerates. 

__Status__: Incomplete

So far, degenerates are handled correctly in some cases, but it is not yet consistent or stable, so more work is needed. Suggestions or contributions are welcome! 

--------------------------------------------------------

# Efficient Clipping of Arbitrary Polygons using OpenGPL

Based on the paper "Efficient Clipping of Arbitrary Polygons" by Günther Greiner (greiner[at]informatik.uni-erlangen.de) and Kai Hormann (hormann[at]informatik.tu-clausthal.de), ACM Transactions on Graphics 1998;17(2):71-83.

Available at: <http://www.inf.usi.ch/hormann/papers/Greiner.1998.ECO.pdf>


## Motivation

This work was created for educational purposes only, as an implementation in Python of the above algorithm, for a class in Graphical Computation.

To study the algorithm, inspect the file `polygon.py`. It can be imported and used in other contexts (i.e., not OpenGL).

### Import

```python
> import polygon
> help(polygon)
> from polygon import *
```


## Command line

The command line interface (`polyclip.py`) is provided for demonstration or testing purposes, using OpenGL.

### Requirements

Supports Python 2.5 or later.

Requires **PyOpenGL** (version 3 as of this writing). If you have pip, install is easy:

`pip install pyopengl`


## Usage

Supported operations are: union, intersection and difference.

### Polygon overrides

Subject and clip polygon can be defined per command line option. Defaults for the subject and clip polygon are set at the beggining of the file for easy edit, but they can be overriden from the command line using the options `--subj-poly` and `--clip-poly`.

**Example:**

`polyclip.py --subj-poly="1.5, 1.25; 7.5, 2.5; 4, 3; 4.5, 6.5"`

### Options

Type `polyclip.py -h` for available options. Press `Esc` to exit.
