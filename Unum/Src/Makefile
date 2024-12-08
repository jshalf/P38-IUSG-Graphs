CXX = g++
 #-DHAVE_GL
PATHBASE = /opt/homebrew
 # /usr/local
CFLAGS = -I$(PATHBASE)/include -I./SoftPosit/source/include 
LINK = -v -L$(PATHBASE)/lib -L. -L./SoftPosit/build/Linux-x86_64-GCC -lgmpxx -lgmp -lsoftposit -lm 
FL = /usr/local
ARCHIVE = /usr/bin/ar cr

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
VPATH = MatrixAlgebra/Sanjana 
# OBJECTS = Posit32.o PositBase.o bfloat16.o gmp_helpers.o 

all: softposit Test2 matrix 

matrix: 
	cd MatrixAlgebra && $(MAKE)

libPosit.a: PositBase.o Posit32.o 
	$(ARCHIVE) libPosit.a PositBase.o Posit32.o 

softposit: 
	cd SoftPosit/build/Linux-x86_64-GCC && $(MAKE)

Test.o: Test.cc PositBase.hh Posit32.hh gmp_helpers.hh
	$(CXX) $(CFLAGS) -c Test.cc

Test: Test.o libPosit.a
	$(CXX) $(CFLAGS) -o Test Test.o $(LINK) -lPosit

Test2: Test2.o SoftPosit $(OBJECTS)
	$(CXX) $(CFLAGS) -o Test2 Test2.o $(OBJECTS) $(LINK) 

#Viewer: Viewer.o ViewerMain.o
#	$(CXX) $(CFLAGS) -o Viewer ${INC} Viewer.o ViewerMain.o ${LINK_FL}

Test2.o: Test2.cc 
	$(CXX) $(CFLAGS) -c Test2.cc

#ViewerMain.o: ViewerMain.cxx Viewer.h
#	$(CXX) $(CFLAGS) -c ${INC} ViewerMain.cxx

#Viewer.o: Viewer.cxx Viewer.h
#	$(CXX) $(CFLAGS) -c ${INC} Viewer.cxx

PositBase.o: PositBase.cc PositBase.hh 
	$(CXX) $(CFLAGS) -c PositBase.cc

Posit32.o: Posit32.cc Posit32.hh 
	$(CXX) $(CFLAGS) -c Posit32.cc 

#Posit32nc.o: Posit32nc.cc Posit32nc.hh 
#	$(CXX) $(CFLAGS) -c Posit32nc.cc

bfloat16.o: bfloat16.cc bfloat16.hh 
	$(CXX) $(CFLAGS) -c bfloat16.cc

gmp_helpers.o: gmp_helpers.cc gmp_helpers.hh
	$(CXX) $(CFLAGS) -c gmp_helpers.cc

#Posit.o: Posit.cc Posit.hh
#	$(CXX) $(CFLAGS) -c Posit.cc

clean:
	rm -f *.o 
	cd MatrixAlgebra && rm -f *.o

distclean: clean
	rm Test Test2 libPosit.a 

