
CXXFLAGS = -std=c++11 -O3 -Wall
LIBS = -lgdal -lshp -lm
OBJS1 = src/screening.o src/model2D.o src/TIFF_IO.o src/reservoir.o src/coordinates.o src/phes_base.o 
OBJS2 = src/pairing.o src/model2D.o src/TIFF_IO.o src/reservoir.o src/coordinates.o src/phes_base.o 
OBJS3 = src/pretty_set.o src/reservoir.o src/model2D.o src/TIFF_IO.o src/coordinates.o src/phes_base.o
OBJS4 = src/constructor.o src/reservoir.o src/model2D.o src/TIFF_IO.o src/coordinates.o src/phes_base.o src/kml.o
INCDIRS = -Iinclude


utils: bin/screening bin/pairing bin/pretty_set bin/constructor

bin/screening: $(OBJS1)
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJS1) $(LIBS) -o $@

bin/pairing: $(OBJS2)
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJS2) $(LIBS) -o $@

bin/pretty_set: $(OBJS3)
	g++  $(CXXFLAGS) $(LDFLAGS) $(OBJS3) $(LIBS) -o $@

bin/constructor: $(OBJS4)
	g++  $(CXXFLAGS) $(LDFLAGS) $(OBJS4) $(LIBS) -o $@

.c.o:
	g++ $(CXXFLAGS) $(INCDIRS) -c -o $@ $<

clean:
	rm -f $(OBJS1) $(OBJS2) bin/screening bin/pairing bin/pretty_set bin/constructor