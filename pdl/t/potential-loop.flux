##############################
# Potential loop system, 3x3 fluxon grid.
#

GLOBAL PHOTOSPHERE 0 0 0 0 0 1

FRAME 1
NEW 1 0 -0 0.001 1
NEW 2 0 1 0.001 -1
LINE 1 1 2 1
VERTEX 1 1 1 0 0.1 0.1
VERTEX 1 2 2 0 0.2 0.1
VERTEX 1 3 3 0 0.3 0.1
VERTEX 1 4 4 0 0.4 0.1
VERTEX 1 5 5 0 0.5 0.1
VERTEX 1 6 6 0 0.6 0.1
VERTEX 1 7 7 0 0.7 0.1
VERTEX 1 8 8 0 0.8 0.1
VERTEX 1 9 9 0 0.9 0.1

NEW 3 0.02 -0 0.001 1
NEW 4 0.02 1 0.001 -1
LINE 2 3 4 1
VERTEX 2 101 1 0.02 0.1 0.1
VERTEX 2 102 2 0.02 0.2 0.1
VERTEX 2 103 3 0.02 0.3 0.1
VERTEX 2 104 4 0.02 0.4 0.1
VERTEX 2 105 5 0.02 0.5 0.1
VERTEX 2 106 6 0.02 0.6 0.1
VERTEX 2 107 7 0.02 0.7 0.1
VERTEX 2 108 8 0.02 0.8 0.1
VERTEX 2 109 9 0.02 0.9 0.1

NEW 5 0.04 -0 0.001 1
NEW 6 0.04 1 0.001 -1
LINE 3 5 6 1
VERTEX 3 201 1 0.04 0.1 0.1
VERTEX 3 202 2 0.04 0.2 0.1
VERTEX 3 203 3 0.04 0.3 0.1
VERTEX 3 204 4 0.04 0.4 0.1
VERTEX 3 205 5 0.04 0.5 0.1
VERTEX 3 206 6 0.04 0.6 0.1
VERTEX 3 207 7 0.04 0.7 0.1
VERTEX 3 208 8 0.04 0.8 0.1
VERTEX 3 209 9 0.04 0.9 0.1

NEW 7 0 -0.02 0.001 1
NEW 8 0 1.02 0.001 -1
LINE 4 7 8 1
VERTEX 4 1001 1 0 0.1 0.2
VERTEX 4 1002 2 0 0.2 0.2
VERTEX 4 1003 3 0 0.3 0.2
VERTEX 4 1004 4 0 0.4 0.2
VERTEX 4 1005 5 0 0.5 0.2
VERTEX 4 1006 6 0 0.6 0.2
VERTEX 4 1007 7 0 0.7 0.2
VERTEX 4 1008 8 0 0.8 0.2
VERTEX 4 1009 9 0 0.9 0.2

NEW 9 0.02 -0.02 0.001 1
NEW 10 0.02 1.02 0.001 -1
LINE 5 9 10 1
VERTEX 5 1101 1 0.02 0.1 0.2
VERTEX 5 1102 2 0.02 0.2 0.2
VERTEX 5 1103 3 0.02 0.3 0.2
VERTEX 5 1104 4 0.02 0.4 0.2
VERTEX 5 1105 5 0.02 0.5 0.2
VERTEX 5 1106 6 0.02 0.6 0.2
VERTEX 5 1107 7 0.02 0.7 0.2
VERTEX 5 1108 8 0.02 0.8 0.2
VERTEX 5 1109 9 0.02 0.9 0.2

NEW 11 0.04 -0.02 0.001 1
NEW 12 0.04 1.02 0.001 -1
LINE 6 11 12 1
VERTEX 6 1201 1 0.04 0.1 0.2
VERTEX 6 1202 2 0.04 0.2 0.2
VERTEX 6 1203 3 0.04 0.3 0.2
VERTEX 6 1204 4 0.04 0.4 0.2
VERTEX 6 1205 5 0.04 0.5 0.2
VERTEX 6 1206 6 0.04 0.6 0.2
VERTEX 6 1207 7 0.04 0.7 0.2
VERTEX 6 1208 8 0.04 0.8 0.2
VERTEX 6 1209 9 0.04 0.9 0.2

NEW 13 0 -0.04 0.001 1
NEW 14 0 1.04 0.001 -1
LINE 7 13 14 1
VERTEX 7 2001 1 0 0.1 0.3
VERTEX 7 2002 2 0 0.2 0.3
VERTEX 7 2003 3 0 0.3 0.3
VERTEX 7 2004 4 0 0.4 0.3
VERTEX 7 2005 5 0 0.5 0.3
VERTEX 7 2006 6 0 0.6 0.3
VERTEX 7 2007 7 0 0.7 0.3
VERTEX 7 2008 8 0 0.8 0.3
VERTEX 7 2009 9 0 0.9 0.3

NEW 15 0.02 -0.04 0.001 1
NEW 16 0.02 1.04 0.001 -1
LINE 8 15 16 1
VERTEX 8 2101 1 0.02 0.1 0.3
VERTEX 8 2102 2 0.02 0.2 0.3
VERTEX 8 2103 3 0.02 0.3 0.3
VERTEX 8 2104 4 0.02 0.4 0.3
VERTEX 8 2105 5 0.02 0.5 0.3
VERTEX 8 2106 6 0.02 0.6 0.3
VERTEX 8 2107 7 0.02 0.7 0.3
VERTEX 8 2108 8 0.02 0.8 0.3
VERTEX 8 2109 9 0.02 0.9 0.3

NEW 17 0.04 -0.04 0.001 1
NEW 18 0.04 1.04 0.001 -1
LINE 9 17 18 1
VERTEX 9 2201 1 0.04 0.1 0.3
VERTEX 9 2202 2 0.04 0.2 0.3
VERTEX 9 2203 3 0.04 0.3 0.3
VERTEX 9 2204 4 0.04 0.4 0.3
VERTEX 9 2205 5 0.04 0.5 0.3
VERTEX 9 2206 6 0.04 0.6 0.3
VERTEX 9 2207 7 0.04 0.7 0.3
VERTEX 9 2208 8 0.04 0.8 0.3
VERTEX 9 2209 9 0.04 0.9 0.3

FINISHED
