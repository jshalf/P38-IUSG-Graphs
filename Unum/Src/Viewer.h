// generated by Fast Light User Interface Designer (fluid) version 1.0304

#ifndef Viewer_h
#define Viewer_h
#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Scroll.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Counter.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Value_Input.H>

class PositViewer {
public:
  PositViewer();
  Fl_Double_Window *mainWindow;
  static Fl_Menu_Item menu_X[];
  static Fl_Menu_Item menu_Y[];
  static Fl_Menu_Item menu_Operator[];
  void show();
};
#endif
