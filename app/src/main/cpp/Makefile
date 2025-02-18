export CC_USE_MATH = 1
export CC_USE_JSMN = 1
export VKK_USE_UI  = 1

TARGET   = atmosphere-elek
CLASSES  = atmo_renderer atmo_uiScreen atmo_uiWindowHud atmo_uiInfoPanel
SOURCE   = $(TARGET).c $(CLASSES:%=%.c)
OBJECTS  = $(TARGET).o $(CLASSES:%=%.o)
HFILES   = $(CLASSES:%=%.h)
OPT      = -O2 -Wall
CFLAGS   = $(OPT) -I. -I$(VULKAN_SDK)/include `sdl2-config --cflags` -DA3D_GL2
LDFLAGS  = -Llibvkk -lvkk -L$(VULKAN_SDK)/lib -lvulkan `sdl2-config --libs` -Llibbfs -lbfs -Llibsqlite3 -lsqlite3 -Llibxmlstream -lxmlstream -Ltexgz -ltexgz -Llibgltf -lgltf -Llibcc -lcc -Llibexpat/expat/lib -lexpat -lm -ldl -lpthread -lz -ljpeg
CCC      = gcc

all: $(TARGET) libcc libgltf libbfs libsqlite3 libxmlstream texgz libvkk libexpat

$(TARGET): $(OBJECTS) libcc libgltf libbfs libsqlite3 libxmlstream texgz libvkk libexpat
	$(CCC) $(OPT) $(OBJECTS) -o $@ $(LDFLAGS)

.PHONY: libcc libgltf libbfs libsqlite3 libxmlstream texgz libvkk libexpat

libcc:
	$(MAKE) -C libcc

libgltf:
	$(MAKE) -C libgltf

libbfs:
	$(MAKE) -C libbfs

libsqlite3:
	$(MAKE) -C libsqlite3

libxmlstream:
	$(MAKE) -C libxmlstream

texgz:
	$(MAKE) -C texgz

libvkk:
	$(MAKE) -C libvkk

libexpat:
	$(MAKE) -C libexpat/expat/lib

clean:
	rm -f $(OBJECTS) *~ \#*\# $(TARGET)
	$(MAKE) -C libcc clean
	$(MAKE) -C libgltf clean
	$(MAKE) -C libbfs clean
	$(MAKE) -C libsqlite3 clean
	$(MAKE) -C libxmlstream clean
	$(MAKE) -C texgz clean
	$(MAKE) -C libvkk clean
	$(MAKE) -C libexpat/expat/lib clean

$(OBJECTS): $(HFILES)
