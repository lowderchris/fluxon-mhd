##############################
# simple_two_fluxons --
# a simple test case with two fluxons of two elements each

FRAME 1
NEW 1	-0.5	-1	00	1
NEW 2	-0.5	1	00	-1

NEW 3   0.5	-1	00      1
NEW 4	0.5	1	00	-1

LINE 1	1 2 1.0
VERTEX 1 11  -0.5 0 0 

LINE 2 3 4 1.0
VERTEX 2 21  0.5 0 0

