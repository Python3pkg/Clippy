# Clippy

Polygon clipping in pure Python. Various implementations including Kim-Kim's degeneracy extension of Greiner-Horman, an abondoned version of the Forster-Overfelt extension of G&H (after discovering their proposed algorithm was flawed to begin with), Vatti, and Martinez et al. Some of these are copied directly from existing libs. Partly for educational purposes and partly for portable pure-Python clipping. 

Karim Bahgat 2016

## Status
No high-level user API yet. Implementations are roughly complete, but seem to contain mistakes as none of the algorithms work correctly. Lots of debugging needed. 