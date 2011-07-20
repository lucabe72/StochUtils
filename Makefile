CFLAGS = -Wall -Wextra -g
CFLAGS += -O8
LDLIBS = -lm

all: meschac avg pmftest pmf-y pmf-yt pmf-rnd matrix pmfgen stoch pmf-iz-transform

avg: avg.o pmf.o pmf-file.o gamma.o y.o
avg-c1: avg-c1.o pmf.o pmf-file.o
pmftest: pmftest.o pmf.o pmf-file.o
pmf-y: pmf-y.o pmf.o pmf-file.o gamma.o y.o v.o pmf-modify.o dl.o
pmf-yt: pmf-yt.o pmf.o pmf-file.o gamma.o z.o y.o v.o pmf-modify.o dl.o
pmf-rnd: pmf-rnd.o pmf.o pmf-file.o pmf-sample.o
matrix: matrix.o pmf.o pmf-file.o cdf.o qbdm.o meschac/meschach.a
pmfgen: pmfgen.o pmf.o pmf-sample.o pmf-modify.o
ifdef GSL
pmfgen: LDFLAGS += `pkg-config --libs gsl` 
pmfgen.o: CPPFLAGS += -DGSL
endif
stoch: stoch.o pmf.o pmf-file.o driver.o pseudo.o generic.o
pmf-iz-transform: pmf-iz-transform.o pmf.o pmf-file.o generic.o
pmf2cdf: pmf2cdf.o pmf.o pmf-file.o cdf.o

meschac/meschach.a:
	make -C meschac
clean:
	rm -f *.o avg avg-c1 pmftest pmf-y pmf-yt pmf-rnd matrix stoch pmf-iz-transform meschach.a
