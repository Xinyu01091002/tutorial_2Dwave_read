# CFLAGS += -O2 -disable-dimensions -grid=quadtree -D_MPI=0
CFLAGS += -O -Og -disable-dimensions -grid=quadtree -D_MPI=0
# CC='mpicc -D_MPI=8' make wave.tst
include $(BASILISK)/Makefile.defs
