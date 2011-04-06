CPPFLAGS = -I..
CFLAGS += -Wall -Wextra -g -O8
ifdef DEBUG
LDFLAGS += -lefence
endif
LDLIBS = -lm

PROGRAM = stoch
OBJS = ../pmf.o ../pmf-file.o driver.o pseudo.o generic.o

all: $(PROGRAM)

stoch: stoch.o $(OBJS)

clean:
	rm -f stoch
	rm -f *.o
