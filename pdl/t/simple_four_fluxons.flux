

######################################

FRAME 1
GLOBAL FORCES b_eqa f_curv_hm f_p_eqa_radial f_vertex
GLOBAL STATE 1  (LOADING)
GLOBAL PHOTOSPHERE <NONE>

NEW	1	-0.5	-1	0	1
NEW	2	-0.5	1	0	-1
NEW	3	0.5	-1	0	1
NEW	4	0.5	1	0	-1
NEW	5	0	-1	-0.5	1
NEW	6	0	1	-0.5	-1
NEW	7	0	-1	0.5	1
NEW	8	0	1	0.5	-1

# FC 1
LINE	1	1	2	1
	VERTEX	1	11	1	-0.5	0	0

# FC 3
LINE	2	3	4	1
	VERTEX	2	21	1	0.5	0	0

# FC 5
LINE	3	5	6	1
	VERTEX	3	31	1	0	0	-0.5

# FC 7
LINE	4	7	8	1
	VERTEX	4	41	1	0	0	0.5

