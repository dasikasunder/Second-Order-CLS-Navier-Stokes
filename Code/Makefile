CC=gcc
CFLAGS=-O3 -Wall -Wextra -MD
SRC := $(wildcard *.c)
LIBS = -lgsl -lm -llis -lopenblas

OBJ=${SRC:.c=.o}

%.o: %.c
	${CC} ${CFLAGS} -o  $@  -c $< ${LIBS}
a.out: ${OBJ}
	${CC} ${CFLAGS} -o $@  ${OBJ} ${LIBS}
clean:
	@rm -f *.o *.d a.out
distclean:
	@rm -f *.o *.d *.csv a.out
scrub:
	@rm -f *.csv

-include $(OBJ:.o=.d)
