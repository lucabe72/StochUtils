CPPFLAGS += -Iinclude
CFLAGS += -Wall -Wextra -g
CFLAGS += -O8
LDLIBS = -lm

VPATH=src:Tools

all: avg avg-c1 confint pmf-y pmf-yt 

avg: avg.o pmf.o pmf-file.o gamma.o y.o
avg-c1: avg-c1.o pmf.o pmf-file.o
pmf-y:  pmf-y.o  pmf.o pmf-file.o gamma.o y.o v.o pmf-modify.o dl.o
pmf-yt: pmf-yt.o pmf.o pmf-file.o gamma.o z.o y.o v.o pmf-modify.o dl.o

clean:
	rm -f *.o avg avg-c1 confint pmf-y pmf-yt
