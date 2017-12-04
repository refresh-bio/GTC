all: gtc

ifdef MSVC     # Avoid the MingW/Cygwin sections
    uname_S := Windows
else                          # If uname not available => 'not' 
    uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
endif

# default install location (binary placed in the /bin folder)
prefix      = /usr/local

# optional install location
exec_prefix = $(prefix)

LIBS_DIR=lib
INCLUDES_DIR=include
CFLAGS=-Wall -O3 -m64 -std=c++11 -pthread -I $(INCLUDES_DIR) -mpopcnt 
CLINK=-O3 -lm -lz -std=c++11 -lpthread -mpopcnt 


ifeq ($(uname_S),Linux)
    CC=g++      
    # check if CPU supports SSE4.2 
    HAVE_SSE4=$(filter-out 0,$(shell grep sse4.2 /proc/cpuinfo | wc -l))
endif
ifeq ($(uname_S),Darwin)
    CC=clang++
    # check if CPU supports SSE4.2
    HAVE_SSE4=$(filter-out 0,$(shell  sysctl machdep.cpu.features| grep SSE4.2 - | wc -l))
endif

CFLAGS+=$(if $(HAVE_SSE4),-msse4.2)
ifeq ($(HAVE_SSE4),)
    CFLAGS+=-msse2
endif

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

gtc:	src/bit_memory.o \
	src/block_init_compressor.o \
	src/buffered_bm.o \
	src/compressed_pack.o \
	src/compression_settings.o \
	src/decompressor.o \
	src/end_compressor.o \
	src/huffman.o \
	src/main.o \
	src/my_vcf.o \
	src/samples.o \
	src/VCFManager.o \
	include/cpp-mmf/memory_mapped_file.o
	$(CC) $(CLINK) -o gtc \
	src/bit_memory.o \
	src/block_init_compressor.o \
	src/buffered_bm.o \
	src/compressed_pack.o \
	src/compression_settings.o \
	src/decompressor.o \
	src/end_compressor.o \
	src/huffman.o \
	src/main.o \
	src/my_vcf.o \
	src/samples.o \
	src/VCFManager.o \
	include/cpp-mmf/memory_mapped_file.o \
	$(LIBS_DIR)/libhts.a \
	$(LIBS_DIR)/libsdsl.a

clean:
	-rm include/cpp-mmf/*.o
	-rm src/*.o
	-rm gtc
	
install:
	mkdir -p -m 755 $(exec_prefix)/bin
	cp gtc $(exec_prefix)/bin/
	
uninstall:
	rm  $(exec_prefix)/bin/gtc
	
