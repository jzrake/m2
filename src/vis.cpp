extern "C" {
#include "m2.h" 
}
#if (M2_HAVE_GLUT)



#include <string>
#include <map>
#include <cmath>
#include "glui/glui.h"


static GLUI *glui;
static int last_x, last_y;
static int main_window;
static float xy_aspect;




class SimulationController
{
private:
  int user_id;
  enum {
    ACTION_RESET_VIEW,
    ACTION_STEP,
    ACTION_QUIT,
  } ;
public:
  float rotationX;
  float rotationY;

  /* live variables */
  float ZoomLevel;
  int AutoAdvance;

  SimulationController(GLUI_Panel *parent);
  void keyboard(int key, int x, int y) { }
  void draw();
  void refresh();
  static std::map<int, SimulationController*> instances;
  static void refresh_cb(int user_id);
  static void action_cb(int action_id);
} ;
std::map<int, SimulationController*> SimulationController::instances;


class DatasetController
{
private:
  double DataRange[2];
  GLfloat *ColorData;
  GLfloat *VertexData;
  double AngularOffset;
  int user_id;
  std::vector<int> obj_keys;

  GLUI_StaticText *text_min;
  GLUI_StaticText *text_max;

  /* live variables */
  int DrawMesh;
  int LogScale;
  int ColorMap;
  int DataMember;
public:
  DatasetController(GLUI_Panel *parent);
  void SetAngularOffset(double x) { AngularOffset = x; }
  void keyboard(int key, int x, int y) { }
  void draw();
  void refresh();
  static std::map<int, DatasetController*> instances;
  static void refresh_cb(int user_id);
} ;
std::map<int, DatasetController*> DatasetController::instances;



static SimulationController *sim_controller = NULL;
static m2sim *M2;
static void color_map(double val, GLfloat rgb[3], int ct);










void myGlutKeyboard(unsigned char Key, int x, int y)
{
  if (Key == 27) {
    exit(0);
  }
  else {
    //    dataset->keyboard(Key, x, y);
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

  glui->sync_live();

  if (sim_controller->AutoAdvance) {
    m2sim_drive(M2);
    sim_controller->refresh();
  }
}


void myGlutMouse(int button, int button_state, int x, int y)
{
  if (button == GLUT_LEFT_BUTTON && button_state == GLUT_DOWN) {
    last_x = x;
    last_y = y;
  }
}
void myGlutMotion(int x, int y)
{
  sim_controller->rotationX += (float) (y - last_y);
  sim_controller->rotationY += (float) (x - last_x);
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
void myGlutDisplay(void)
{
  sim_controller->draw();
  glutSwapBuffers(); 
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


  /* Build the GUI */
  glui = GLUI_Master.create_glui("GLUI", 0, 768+50, 50);
  new GLUI_StaticText(glui, "M2 controls"); 
  GLUI_Panel *dataset_panel = new GLUI_Panel(glui, "data sets");
  DatasetController *ds1 = new DatasetController(dataset_panel);
  new GLUI_Column(dataset_panel);
  DatasetController *ds2 = new DatasetController(dataset_panel);
  GLUI_Panel *general_panel = new GLUI_Panel(glui, "simulation controls");
  sim_controller = new SimulationController(general_panel);

  ds1->SetAngularOffset(M2_PI);
  ds2->SetAngularOffset(0.0);

  sim_controller->refresh();

  glui->set_main_gfx_window(main_window);
  GLUI_Master.set_glutIdleFunc(myGlutIdle);
  glutMainLoop();
}




SimulationController::SimulationController(GLUI_Panel *parent)
  : ZoomLevel(0.0),
    AutoAdvance(0)
{
  user_id = instances.size();
  instances[user_id] = this;
  new GLUI_Checkbox(parent, "Auto advance", &AutoAdvance);
  GLUI_Scrollbar *zoom_sb = new GLUI_Scrollbar(parent,
					       "zoom level",
					       GLUI_SCROLL_HORIZONTAL,
					       &ZoomLevel);
  new GLUI_Button(parent, "Reset view", ACTION_RESET_VIEW, action_cb);
  new GLUI_Button(parent, "Step", ACTION_STEP, action_cb);
  new GLUI_Button(parent, "Quit", ACTION_QUIT, action_cb);
  zoom_sb->set_float_limits(-2.0, 2.0);
}
void SimulationController::refresh_cb(int user_id)
{
  instances[user_id]->refresh();
}
void SimulationController::action_cb(int action_id)
{
  switch (action_id) {
  case ACTION_RESET_VIEW:
    sim_controller->ZoomLevel = 0.0;
    sim_controller->rotationX = 0.0;
    sim_controller->rotationY = 0.0;
    break;
  case ACTION_STEP:
    m2sim_drive(M2);
    sim_controller->refresh();
    break;
  case ACTION_QUIT:
    exit(0);
    break;
  }
}
void SimulationController::refresh()
{
  std::map<int, DatasetController*>::iterator d =
    DatasetController::instances.begin();
  while (d != DatasetController::instances.end()) {
    d->second->refresh();
    ++d;
  }
}
void SimulationController::draw()
{
  GLfloat zoom = powf(10.0, ZoomLevel);

  glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
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


  std::map<int, DatasetController*>::iterator d =
    DatasetController::instances.begin();
  while (d != DatasetController::instances.end()) {
    d->second->draw();
    ++d;
  }


  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0.0, 100.0, 0.0, 100.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glColor3ub(0, 0, 0);
  glRasterPos2i(10, 10);

  std::string text = "test";

  for (unsigned int i=0; i<text.length(); ++i) {
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, text[i]);
  }
}


DatasetController::DatasetController(GLUI_Panel *parent)
  :
  ColorData(NULL),
  VertexData(NULL),
  AngularOffset(0.0),

  /* live variables */
  DrawMesh(0),
  LogScale(0),
  ColorMap(3),
  DataMember(0)
{
  user_id = instances.size();
  instances[user_id] = this;
  GLUI_RadioGroup *radio = new GLUI_RadioGroup(parent, &DataMember, user_id, refresh_cb);
  new GLUI_RadioButton(radio, "mass density"); obj_keys.push_back(M2_COMOVING_MASS_DENSITY);
  new GLUI_RadioButton(radio, "gas pressure"); obj_keys.push_back(M2_GAS_PRESSURE);
  new GLUI_RadioButton(radio, "v1"); obj_keys.push_back(M2_VELOCITY1);
  new GLUI_RadioButton(radio, "v2"); obj_keys.push_back(M2_VELOCITY2);
  new GLUI_RadioButton(radio, "v3"); obj_keys.push_back(M2_VELOCITY3);
  new GLUI_RadioButton(radio, "B1"); obj_keys.push_back(M2_MAGNETIC1);
  new GLUI_RadioButton(radio, "B2"); obj_keys.push_back(M2_MAGNETIC2);
  new GLUI_RadioButton(radio, "B3"); obj_keys.push_back(M2_MAGNETIC3);
  new GLUI_Separator(parent);
  text_min = new GLUI_StaticText(parent, "min");
  text_max = new GLUI_StaticText(parent, "max");
  new GLUI_Separator(parent);
  new GLUI_Checkbox(parent, "draw mesh", &DrawMesh, user_id, refresh_cb);
  new GLUI_Checkbox(parent, "log scale", &LogScale, user_id, refresh_cb);
  GLUI_Spinner *cm = new GLUI_Spinner (parent, "color map", &ColorMap, user_id, refresh_cb);
  cm->set_int_limits(0, 5);
  this->refresh();
}
void DatasetController::draw()
{
  for (int n=0; n<M2->local_grid_size[0]; ++n) {
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
}
void DatasetController::refresh()
{
  int *L = M2->local_grid_size;
  double y;
  double x00[4]; /* vertex coordinates */
  double x01[4];
  double x11[4];
  double x10[4];
  double x00c[4]; /* in cartesian geometry */
  double x01c[4];
  double x11c[4];
  double x10c[4];
  m2vol *V;
  VertexData = (GLfloat*) realloc(VertexData, L[0]*3*4*sizeof(GLfloat));
  ColorData = (GLfloat*) realloc(ColorData, L[0]*3*sizeof(GLfloat));
  for (int n=0; n<L[0]; ++n) {
    V = M2->volumes + n;
    y = m2aux_get(&V->aux, obj_keys[DataMember]);
    if (LogScale) {
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

  char min[256];
  char max[256];
  snprintf(min, 256, "min: %3.2f", DataRange[0]);
  snprintf(max, 256, "max: %3.2f", DataRange[1]);
  text_min->set_text(min);
  text_max->set_text(max);

  for (int n=0; n<L[0]; ++n) {
    V = M2->volumes + n;
    x00[1] = V->x0[1];
    x00[2] = V->x0[2];
    x00[3] = V->x0[3] + AngularOffset;
    x01[1] = V->x1[1];
    x01[2] = V->x0[2];
    x01[3] = V->x0[3] + AngularOffset;
    x11[1] = V->x1[1];
    x11[2] = V->x1[2];
    x11[3] = V->x0[3] + AngularOffset;
    x10[1] = V->x0[1];
    x10[2] = V->x1[2];
    x10[3] = V->x0[3] + AngularOffset;
    m2_to_cartesian(x00, x00c, M2->geometry);
    m2_to_cartesian(x01, x01c, M2->geometry);
    m2_to_cartesian(x11, x11c, M2->geometry);
    m2_to_cartesian(x10, x10c, M2->geometry);
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
    y = m2aux_get(&V->aux, obj_keys[DataMember]);
    if (LogScale) {
      y = log10(fabs(y));
    }
    y -= DataRange[0];
    y /= DataRange[1] - DataRange[0];
    color_map(y, &ColorData[3*n], ColorMap);
  }
}
void DatasetController::refresh_cb(int user_id)
{
  instances[user_id]->refresh();
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
