extern "C" {
#include "m2.h" 
}
#if (M2_HAVE_GLUT)

/****************************************************************************

  example3.cpp

  A GLUT program using all the features of the GLUI User Interface Library
  (except columns, featured in example4.cpp)

  -----------------------------------------------------------------------
	   
  9/9/98 Paul Rademacher (rademach@cs.unc.edu)

****************************************************************************/

#include <string>
#include <map>
#include <cmath>
#include "glui/glui.h"

int   last_x, last_y;
float xy_aspect;
float rotationX = 0.0, rotationY = 0.0;


/** These are the live variables passed into GLUI ***/
int   obj_type = 0;
int   ColorTable = 3;
int   main_window;
int   counter = 0;
int   AutoPlay = 0;
int   LogScaling = 0;
int   DrawMesh = 0;
float scale = 0.0;
std::string text = "";


/********** User IDs for callbacks ********/
#define OPEN_CONSOLE_ID      100
#define CMD_HIST_RESET_ID    101
#define CMD_CLOSE_ID         102
#define LIGHT0_ENABLED_ID    200
#define LIGHT1_ENABLED_ID    201
#define LIGHT0_INTENSITY_ID  250
#define LIGHT1_INTENSITY_ID  251
#define CHANGE_COLORTABLE_ID 350
#define CHANGE_FIELD_ID      351
#define CHANGE_LOGSCALING_ID 352
#define CHANGE_DRAWMESH_ID   353



static GLUI *glui;
class M2GL_Dataset;
static m2sim *M2;
static std::vector<int> obj_keys;
static void color_map(double val, GLfloat rgb[3], int ct);
static M2GL_Dataset *dataset;


class M2GL_Dataset
{
private:
  m2sim *m2;
  int DataMember;
  //  int ColorTable;
  //  int DrawMesh;
  //  int LogScaling;
  double DataRange[2];
  GLfloat *VertexData;
  GLfloat *ColorData;
public:
  M2GL_Dataset(m2sim *m2) :
    m2(m2),
    DataMember(M2_COMOVING_MASS_DENSITY),
    //    ColorTable(3),
    //    DrawMesh(0),
    //    LogScaling(0),
    VertexData(NULL),
    ColorData(NULL)
  { }
  void reload_data()
  {
    int n;
    int *L = m2->local_grid_size;
    m2vol *V;
    double y;
    double x00[4]; /* vertex coordinates */
    double x01[4];
    double x11[4];
    double x10[4];
    double x00c[4]; /* in cartesian geometry */
    double x01c[4];
    double x11c[4];
    double x10c[4];


    VertexData = (GLfloat*) realloc(VertexData, L[0]*3*4*sizeof(GLfloat));
    ColorData = (GLfloat*) realloc(ColorData, L[0]*3*sizeof(GLfloat));
    DataMember = M2_COMOVING_MASS_DENSITY;//obj_keys[obj_type];


    for (n=0; n<L[0]; ++n) {
      V = m2->volumes + n;

      y = m2aux_get(&V->aux, DataMember);

      if (LogScaling) {
	if (y < 1e-6) y = 1e-6;
	y = log10(fabs(y));
      }

      if (n == 0) {
	DataRange[0] = y;
	DataRange[1] = y;
      }
      else {
	if (y < DataRange[0]) DataRange[0] = y;
	if (y > DataRange[1]) DataRange[1] = y;
      }
    }

    for (n=0; n<L[0]; ++n) {
      V = m2->volumes + n;

      x00[1] = V->x0[1];
      x00[2] = V->x0[2];
      x00[3] = V->x0[3];
      x01[1] = V->x1[1];
      x01[2] = V->x0[2];
      x01[3] = V->x0[3];
      x11[1] = V->x1[1];
      x11[2] = V->x1[2];
      x11[3] = V->x0[3];
      x10[1] = V->x0[1];
      x10[2] = V->x1[2];
      x10[3] = V->x0[3];

      m2_to_cartesian(x00, x00c, m2->geometry);
      m2_to_cartesian(x01, x01c, m2->geometry);
      m2_to_cartesian(x11, x11c, m2->geometry);
      m2_to_cartesian(x10, x10c, m2->geometry);

      VertexData[4*3*n + 0*3 + 0] = x00c[1];
      VertexData[4*3*n + 0*3 + 1] = x00c[2];
      VertexData[4*3*n + 0*3 + 2] = x00c[3];
      VertexData[4*3*n + 1*3 + 0] = x01c[1];
      VertexData[4*3*n + 1*3 + 1] = x01c[2];
      VertexData[4*3*n + 1*3 + 2] = x01c[3];
      VertexData[4*3*n + 2*3 + 0] = x11c[1];
      VertexData[4*3*n + 2*3 + 1] = x11c[2];
      VertexData[4*3*n + 2*3 + 2] = x11c[3];
      VertexData[4*3*n + 3*3 + 0] = x10c[1];
      VertexData[4*3*n + 3*3 + 1] = x10c[2];
      VertexData[4*3*n + 3*3 + 2] = x10c[3];

      y = m2aux_get(&V->aux, DataMember);

      if (LogScaling) {
	y = log10(fabs(y));
      }

      y -= DataRange[0];
      y /= DataRange[1] - DataRange[0];
      color_map(y, &ColorData[3*n], ColorTable);
    }
  }
  void draw()
  {
    glDisable(GL_LIGHTING);
    for (int n=0; n<m2->local_grid_size[0]; ++n) {
      if (DrawMesh) {
	glBegin(GL_LINE_LOOP);
      }
      else {
	glBegin(GL_QUADS);
      }
      glColor3fv(&ColorData[3*n]);
      glVertex3fv(&VertexData[4*3*n + 0*3]);
      glVertex3fv(&VertexData[4*3*n + 1*3]);
      glVertex3fv(&VertexData[4*3*n + 2*3]);
      glVertex3fv(&VertexData[4*3*n + 3*3]);

      glEnd();
    }
    glEnable(GL_LIGHTING);
  }
  void keyboard(int key, int x, int y)
  {
    switch (key) {
    case 'f':
      ++obj_type %= obj_keys.size();
      this->reload_data();
      GLUI_Master.sync_live_all();
      break;
    case 'm':
      DrawMesh ^= 1;
      GLUI_Master.sync_live_all();
      break;
    case 'l':
      LogScaling ^= 1;
      this->reload_data();
      GLUI_Master.sync_live_all();
      break;
    }
  }
} ;
void myGlutKeyboard(unsigned char Key, int x, int y)
{
  if (Key == 27) {
    exit(0);
  }
  else {
    dataset->keyboard(Key, x, y);
  }
  glutPostRedisplay();
}
void myGlutMenu(int value)
{
  myGlutKeyboard(value, 0, 0);
}
void myGlutIdle(void)
{
  /* According to the GLUT specification, the current window is 
     undefined during an idle callback.  So we need to explicitly change
     it if necessary */
  if (glutGetWindow() != main_window) {
    glutSetWindow(main_window);
  }
  glutPostRedisplay();
  /****************************************************************/
  /*            This demonstrates GLUI::sync_live()               */
  /*   We change the value of a variable that is 'live' to some   */
  /*   control.  We then call sync_live, and the control          */
  /*   associated with that variable is automatically updated     */
  /*   with the new value.  This frees the programmer from having */
  /*   to always remember which variables are used by controls -  */
  /*   simply change whatever variables are necessary, then sync  */
  /*   the live ones all at once with a single call to sync_live  */
  /****************************************************************/

  counter++;
  glui->sync_live();

  if (AutoPlay) {
    m2sim_drive(M2);
    dataset->reload_data();
  }
}


void myGlutMouse(int button, int button_state, int x, int y)
{
  if ( button == GLUT_LEFT_BUTTON && button_state == GLUT_DOWN ) {
    last_x = x;
    last_y = y;
  }
}
void myGlutMotion(int x, int y)
{
  rotationX += (float) (y - last_y);
  rotationY += (float) (x - last_x);
  last_x = x;
  last_y = y;
  glutPostRedisplay(); 
}
void myGlutReshape(int x, int y)
{
  xy_aspect = (float) x / (float) y;
  glViewport(0, 0, x, y);
  glutPostRedisplay();
}
void myGlutDisplay( void )
{
  GLfloat zoom = powf(10.0, scale);
  glClearColor(.9f, .9f, .9f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glFrustum(-xy_aspect*0.08, xy_aspect*0.08, -0.08, 0.08, 0.1, 15.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(0.0, 0.0, -1.0f);
  glRotatef(rotationX, 1.0, 0.0, 0.0);
  glRotatef(rotationY, 0.0, 1.0, 0.0);
  glRotated(-90.0, 1.0, 0.0, 0.0);
  glScalef(zoom, zoom, zoom);
 

  dataset->draw();


  /* Disable lighting and set up ortho projection to render text */
  glDisable(GL_LIGHTING);  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0.0, 100.0, 0.0, 100.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glColor3ub(0, 0, 0);
  glRasterPos2i(10, 10);
  for (unsigned int i=0; i<text.length(); ++i) {
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, text[i]);
  }
  glEnable(GL_LIGHTING);

  glutSwapBuffers(); 
}



class MyClass
{
private:
  int user_id;
  int val1, val2, val3, val4;
public:
  MyClass(GLUI_Panel *parent);
  void keyboard(int key, int x, int y) { }
  void draw() { }
  void refresh() { }
  static std::map<int, MyClass*> instances;
  static void refresh_cb(int user_id)
  {
    instances[user_id]->refresh();
  }
} ;
std::map<int, MyClass*> MyClass::instances;
MyClass::MyClass(GLUI_Panel *parent)
{
  user_id = instances.size();
  instances[user_id] = this;
  new GLUI_Column(parent);
  new GLUI_Checkbox(parent, "check1", &val1, user_id, refresh_cb);
  new GLUI_Checkbox(parent, "check2", &val2, user_id, refresh_cb);
  new GLUI_Checkbox(parent, "check3", &val3, user_id, refresh_cb);
  new GLUI_Checkbox(parent, "check4", &val4, user_id, refresh_cb);
}




class SimulationController
{
private:
  int user_id;
  int val1, val2, val3, val4;
public:
  SimulationController(GLUI_Panel *parent);
  void keyboard(int key, int x, int y) { }
  void draw() { }
  void refresh() { }
  static std::map<int, SimulationController*> instances;
  static void refresh_cb(int user_id)
  {
    instances[user_id]->refresh();
  }
} ;
std::map<int, SimulationController*> SimulationController::instances;
SimulationController::SimulationController(GLUI_Panel *parent)
{
  user_id = instances.size();
  instances[user_id] = this;
  new GLUI_Column(parent);
  new GLUI_Checkbox(parent, "check1", &val1, user_id, refresh_cb);
  new GLUI_Checkbox(parent, "check2", &val2, user_id, refresh_cb);
  new GLUI_Checkbox(parent, "check3", &val3, user_id, refresh_cb);
  new GLUI_Checkbox(parent, "check4", &val4, user_id, refresh_cb);
  new GLUI_Button(parent, "Quit", 0, exit);
}








void m2sim_visualize(m2sim *m2, int argc, char **argv)
{
  M2 = m2;
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowPosition(50, 50);
  glutInitWindowSize(768, 768);
 
  main_window = glutCreateWindow("M2");
  glutDisplayFunc(myGlutDisplay);
  glutReshapeFunc(myGlutReshape);  
  glutKeyboardFunc(myGlutKeyboard);
  glutMotionFunc(myGlutMotion);
  glutMouseFunc(myGlutMouse);

  glEnable(GL_DEPTH_TEST);


  glui = GLUI_Master.create_glui("GLUI", 0, 768+50, 50);
  new GLUI_StaticText(glui, "M2 controls"); 


  GLUI_Panel *dataset_panel = new GLUI_Panel(glui, "data sets");
  new MyClass(dataset_panel);
  new MyClass(dataset_panel);


  GLUI_Panel *general_panel = new GLUI_Panel(glui, "simulation controls");
  new SimulationController(general_panel);


  dataset = new M2GL_Dataset(m2);
  dataset->reload_data();



  glui->set_main_gfx_window(main_window);
  GLUI_Master.set_glutIdleFunc(myGlutIdle);
  glutMainLoop();
}


void color_map(double val, GLfloat rgb[3], int ct)
{
  double rrr, ggg, bbb;
  if( ct == 0 ){
    double nexp = 8.0;
    rrr = exp(-nexp*pow(val-5./6.,2.0)) + .25*exp(-nexp*pow(val+1./6.,2.0));
    ggg = exp(-nexp*pow(val-3./6.,2.0));
    bbb = exp(-nexp*pow(val-1./6.,2.0)) + .25*exp(-nexp*pow(val-7./6.,2.0));
  }else if(ct == 1){
    if( val < .1 ){
      bbb = 4.*(val+.15);
      ggg = 0.0;
      rrr = 0.0;
    }else if( val < .35){
      bbb = 1.0;
      ggg = 4.*(val-.1);
      rrr = 0.0;
    }else if( val < .6 ){
      bbb = 4.*(.6-val);
      ggg = 1.;
      rrr = 4.*(val-.35);
    }else if( val < .85){
      bbb = 0.0;
      ggg = 4.*(.85-val);
      rrr = 1.;
    }else{
      bbb = 0.0;
      ggg = 0.0;
      rrr = 4.*(1.1-val);
    }
  }else if(ct == 2){
    rrr = 2.*val;
    ggg = 1.2*val;
    bbb = .8*val;
  }else if(ct == 3){
    double gam = .8;
    double Amp;
    double r0,g0,b0;
    double hi,lo,x1,x2,x3,x4;
    hi = .8;
    lo = .1;
    if( val > hi ) Amp = .3 + .7*(1.-val)/(1.-hi);
    else if( val < lo ) Amp = .3 + .7*(val)/(lo);
    else Amp = 1.0;
    x1 = .5;
    x2 = .325;
    x3 = .15;
    x4 = 0.;
    if( val > x1 )      r0 = 1.;
    else if( val > x2 ) r0 = (val-x2)/(x1-x2);
    else if( val > x3 ) r0 = 0.;
    else if( val > x4 ) r0 = (val-x3)/(x4-x3);
    else                r0 = 1.;
    x1 = .6625;
    x2 = .5;
    x3 = .275;
    x4 = .15;
    if( val > x1 )      g0 = 0.;
    else if( val > x2 ) g0 = (val-x1)/(x2-x1);
    else if( val > x3 ) g0 = 1.;
    else if( val > x4 ) g0 = (val-x4)/(x3-x4);
    else                g0 = 0.;
    x1 = .325;
    x2 = .275;
    if( val > x1 )      b0 = 0.;
    else if( val > x2 ) b0 = (val-x1)/(x2-x1);
    else                b0 = 1.;
    rrr = pow(Amp*r0,gam);
    ggg = pow(Amp*g0,gam);
    bbb = pow(Amp*b0,gam);
  }else if(ct == 4){
    if( val < .1 ){
      bbb = 4.*(val+.125);
      ggg = 0.0;
      rrr = 0.0;
    }else if( val < .375){
      bbb = 1.0;
      ggg = 4.*(val-.125);
      rrr = 0.0;
    }else if( val < .625 ){
      bbb = 4.*(.625-val);
      rrr = 4.*(val-.375);
      ggg = bbb;
      if( rrr > bbb ) ggg = rrr;
    }else if( val < .875){
      bbb = 0.0;
      ggg = 4.*(.875-val);
      rrr = 1.;
    }else{
      bbb = 0.0;
      ggg = 0.0;
      rrr = 4.*(1.125-val);
    }
  }else if(ct == 5){
    rrr = val;
    ggg = val;
    bbb = val;
  }else{
    rrr = 1.0;
    ggg = 1.0;
    bbb = 1.0;
  }
  rgb[0] = rrr;
  rgb[1] = ggg;
  rgb[2] = bbb;
}


#else
void m2sim_visualize(m2sim *m2, int argc, char **argv)
{

}
#endif /* M2_HAVE_GLUT */
