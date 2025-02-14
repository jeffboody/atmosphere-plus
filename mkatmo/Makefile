export CC_USE_MATH = 1
export CC_USE_JSMN = 1

TARGET   = mkatmo
CLASSES  =
SOURCE   = $(TARGET).c $(CLASSES:%=%.c)
OBJECTS  = $(TARGET).o $(CLASSES:%=%.o)
HFILES   = $(CLASSES:%=%.h)
OPT      = -O2 -Wall
CFLAGS   = $(OPT) -I.
LDFLAGS  = -Ltexgz -ltexgz -Llibcc -lcc -lz -lm
CCC      = gcc

all: $(TARGET) libcc texgz

$(TARGET): $(OBJECTS) libcc texgz
	$(CCC) $(OPT) $(OBJECTS) -o $@ $(LDFLAGS)

.PHONY: libcc texgz

libcc:
	$(MAKE) -C libcc

texgz:
	$(MAKE) -C texgz

clean:
	rm -f $(OBJECTS) *~ \#*\# $(TARGET)
	$(MAKE) -C libcc clean
	$(MAKE) -C texgz clean
	rm libcc texgz lodepng

$(OBJECTS): $(HFILES)
