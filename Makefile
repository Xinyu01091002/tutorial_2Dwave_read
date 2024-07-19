#CFLAGS += -O2 -disable-dimensions -grid=quadtree -D_MPI=0
CFLAGS += -O -Og -disable-dimensions -grid=quadtree -D_MPI=0
include $(BASILISK)/Makefile.defs
