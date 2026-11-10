CC = gcc
CFLAGS = -Wall
LIBS = -lm
LFLAGS = 

DEPS = sbls.h
OBJ := $(DEPS:%.h=%.o)
BINS = period_find

all: ${BINS}

%.o: %.c ${DEPS}
	${CC} -c -o $@ $< ${CFLAGS} ${IFLAGS}

period_find: period_find.o ${OBJ}
	${CC} -o $@ $^ ${CFLAGS} ${LFLAGS} ${LIBS}

clean:
	rm -rf *.o period_find
