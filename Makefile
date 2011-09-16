CC = gcc
EXE = run
CFLAGS = -g -O0 -Wall
GSL_LIBS = -L/usr/local/lib -lgsl -lgslcblas
GSL_INC = -I/usr/local/include/gsl
LDFLAGS = -g -O0 -Wall $(GSL_LIBS)

OBJECTS = main.o

all: $(OBJECTS) $(EXE)

$(EXE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(OJBECTS): %.o: %.c
	$(CC) -c $(CFLAGS) $(GSL_INC) $< -o $@

clean:
	rm -rf $(EXE) $(OBJECTS)
