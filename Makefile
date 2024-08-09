# CFLAGS += -O2 -disable-dimensions -grid=quadtree
CFLAGS += -O -Og -disable-dimensions -grid=quadtree 
# CC='mpicc -D_MPI=8' make wave.tst
include $(BASILISK)/Makefile.defs
