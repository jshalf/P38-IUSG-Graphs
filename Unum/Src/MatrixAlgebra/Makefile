CXX = c++ -g -DHAVE_GL
PATHBASE = /opt/homebrew
# was /usr/local
CFLAGS = -I$(PATHBASE)/include -I.. 
LINK = -L$(PATHBASE)/lib -L/usr/X11R6/lib -L. -L.. -L../SoftPosit/build/Linux-x86_64-GCC -lgmpxx \
	-lpthread -lX11 -lgmp -lm -lboost_system -lboost_filesystem
FL = /usr/local

LDLIBS          = -lpthread  -framework Carbon -framework Cocoa -framework ApplicationServices -framework AudioToolbox
GLDLIBS         = -framework AGL -framework OpenGL -lpthread  -framework Carbon -framework Cocoa -framework ApplicationServices -framework AudioToolbox
LINKFLTK        = ${FL}/lib/libfltk.a
LINKFLTKGL      = ${FL}/lib/libfltk_gl.a
LINKFLTKFORMS   = ${FL}/lib/libfltk_forms.a ${FL}/lib/libfltk.a
LINKFLTKIMG     = ${FL}/lib/libfltk_images.a ${FL}/lib/libfltk.a $(IMAGELIBS)
LINKSHARED      =   ${FL}/lib/libfltk_images.a ${FL}/lib/libfltk_forms.a ${FL}/lib/libfltk.a $(IMAGELIBS)
IMAGELIBS       = -L${FL}/lib -lfltk_png -lz -lfltk_jpeg

FLINC = -I${FL} -I${FL}/FL -I${FL}/include/FL -I${FL}/include/Fl
FLLIB = -L${FL}/lib -lfltk -framework Carbon -framework ApplicationServices
GLLIBS = -lfltk_gl -framework AGL -framework OpenGL -framework AGL -framework OpenGL -lpthread  -framework Carbon -framework Cocoa -framework ApplicationServices -framework AudioToolbox
#LINK = -L$(FL)/lib -lfltk_forms -lfltk_images -lfltk_png -lfltk_jpeg -lfltk_gl -lfltk -framework AGL -framework OpenGL  -framework Carbon -framework ApplicationServices -framework Cocoa -framework AudioToolbox

LINK_FL= ${LINKFLTKGL} ${LINKFLTK} ${GLDLIBS}
XLIBS = -L/usr/X11R6/lib -lX11 -lXmu -lXext
# LIBS = ${FLLIB} ${GLLIBS}
INC = ${FLINC}
#APLRES = /Developer/Tools/Rez -t APPL -o
#RESPATH = $(FL)/FL/mac.r
#GLDLIBS = -framework AGL -framework OpenGL -framework Carbon -framework ApplicationServices
VPATH = ..

OBJ = Matrix.o SuperMatrix.o Quire.o FFT.o helpers.o mmio.o \
	../Posit32.o ../PositBase.o ../bfloat16.o ../gmp_helpers.o #../Posit32nc.o

OBJS = Matrix.o SuperMatrix.o Quire.o FFT.o helpers.o mmio.o \
	Posit32.o PositBase.o bfloat16.o gmp_helpers.o #Posit32nc.o 

HEADERS = $(OBJS:.o=.hh)

all: above Benchmarks Test

p38: Main

above:
	cd .. && $(MAKE) softposit Test2

Test: Test.o Matrix.hh SuperMatrix.hh
	$(CXX) $(CFLAGS) -o Test Test.o $(OBJ) $(LINK)

Main: Main.o Matrix.o mmio.o SuperMatrix.o
	$(CXX) $(CFLAGS) -o Main Main.o mmio.o Matrix.o $(LINK)

Main.o: Main.cc Matrix.hh
	$(CXX) $(CFLAGS) -c Main.cc

Benchmarks: Benchmarks.o $(OBJS) 
	$(CXX) $(CFLAGS) -o Benchmarks Benchmarks.o $(OBJ) $(LINK)

Benchmarks.o: Benchmarks.cc $(HEADERS)
	$(CXX) $(CFLAGS) -c Benchmarks.cc

Test.o: Test.cc $(HEADERS)
	$(CXX) $(CFLAGS) -c Test.cc

Matrix.o: Matrix.cc Matrix.hh $(HEADERS)
	$(CXX) $(CFLAGS) -c Matrix.cc

SuperMatrix.o: SuperMatrix.cc Matrix.hh $(HEADERS)
	$(CXX) $(CFLAGS) -c SuperMatrix.cc

mmio.o: mmio.cc mmio.hh 
	$(CXX) $(CFLAGS) -c mmio.cc

Viewer: Viewer.o ViewerMain.o 
	$(CXX) $(CFLAGS) -o Viewer ${INC} Viewer.o ViewerMain.o ${LINK_FL}

ViewerMain.o: ViewerMain.cxx Viewer.h
	$(CXX) $(CFLAGS) -c ${INC} ViewerMain.cxx

Viewer.o: Viewer.cxx Viewer.h
	$(CXX) $(CFLAGS) -c ${INC} Viewer.cxx

Quire.o: Quire.cc Quire.hh $(HEADERS)
	$(CXX) $(CFLAGS) -c Quire.cc

helpers.o: helpers.cc helpers.hh $(HEADERS)
	$(CXX) $(CFLAGS) -c helpers.cc 

#gmp_helpers.o: gmp_helpers.cc gmp_helpers.hh
#	$(CXX) $(CFLAGS) -c gmp_helpers.cc

FFT.o: FFT.cc FFT.hh $(HEADERS)
	$(CXX) $(CFLAGS) -c FFT.cc

clean:
	rm -f *.o 

distclean: clean
	rm Test Benchmarks

