CFLAGS = -Wall -Wextra -g #-O8
LDLIBS = -lm

all: avg pmftest pmf-y pmf-yt pmf-rnd

avg: avg.o pmf.o pmf-file.o gamma.o y.o
avg-c1: avg-c1.o pmf.o pmf-file.o
pmftest: pmftest.o pmf.o pmf-file.o
pmf-y: pmf-y.o pmf.o pmf-file.o gamma.o y.o v.o pmf-modify.o dl.o
pmf-yt: pmf-yt.o pmf.o pmf-file.o gamma.o z.o y.o v.o pmf-modify.o dl.o
pmf-rnd: pmf-rnd.o pmf.o pmf-file.o pmf-sample.o

clean:
	rm -f *.o avg avg-c1 pmftest pmf-y pmf-yt pmf-rnd
