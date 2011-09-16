CC = gcc
EXE = run

GSL_INC = -I/home/friesen/include
GSL_LIBS = -L/home/friesen/lib -lgsl -lgslcblas

CFLAGS = -g -O0 -Wall $(GSL_INC)
LDFLAGS = -g -O0 -Wall $(GSL_LIBS) -lm

OBJECTS = main.o

all: $(OBJECTS) $(EXE)

$(EXE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(OJBECTS): %.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -rf $(EXE) $(OBJECTS) *.out
