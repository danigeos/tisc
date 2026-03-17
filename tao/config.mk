#------------------------------------------------------------------------
VERSION = tAo_2026-03-15_factorized
#You may need to modify these variables
CC	= gcc #gcc cc
LIBS	= -L$(LIB) -lm -lc
#Options depending on the compiler:
OPTS_mac = -g -w #-Wuninitialized
OPTS_linux = ($OPTS_mac)
OPTS_AIX_RS6000 = -g -O3 #-Q -qsrcmsg #-v
OPTS_SUN = -g #-O3 #-xO5 #-O #Sometimes Sun OS has obscure segmentation problems with -O (like those in grav. anom.)
#Choose your system:
OPTS	= $(OPTS_mac)
#------------------------------------------------------------------------


