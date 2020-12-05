#include <stdio.h>
#include <stdlib.h>
#include "Viewer.h"
/* sample main for volmap */
int main(int argc, char **argv) {

    //RawTCPserver server(5558);
    PositViewer *viewer = new PositViewer;
    
    //Initial global objects.
    
    Fl::visual(FL_DOUBLE|FL_INDEX);
    
    viewer->show();
    return Fl::run();

}

