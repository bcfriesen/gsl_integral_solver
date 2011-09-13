CC = gcc
OBJECTS = main.o
EXE = run
CFLAGS = -g -O0 -Wall
GSL_LIBS = -lgsl -lgslcblas
LDFLAGS = -g -O0 -Wall $(GSL_LIBS)

all: $(OBJECTS) $(EXE)

$(EXE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(OJBECTS): %.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -rf $(EXE) $(OBJECTS)
