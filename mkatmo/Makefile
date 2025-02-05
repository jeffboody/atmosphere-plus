export CC_USE_MATH = 1
export CC_USE_JSMN = 1

TARGET   = mkatmo
CLASSES  =
SOURCE   = $(TARGET).c $(CLASSES:%=%.c)
OBJECTS  = $(TARGET).o $(CLASSES:%=%.o)
HFILES   = $(CLASSES:%=%.h)
OPT      = -O2 -Wall
CFLAGS   = $(OPT) -I.
LDFLAGS  = -Llibcc -lcc -lm
CCC      = gcc

all: $(TARGET) libcc

$(TARGET): $(OBJECTS) libcc
	$(CCC) $(OPT) $(OBJECTS) -o $@ $(LDFLAGS)

.PHONY: libcc

libcc:
	$(MAKE) -C libcc

clean:
	rm -f $(OBJECTS) *~ \#*\# $(TARGET)
	$(MAKE) -C libcc clean
	rm libcc

$(OBJECTS): $(HFILES)
