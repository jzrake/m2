extern "C" {
#include "m2.h" 
}
#if (M2_HAVE_GLUT)


#include <string>
#include <map>
#include <cmath>
#include "glui/glui.h"


static GLUI *glui=NULL;
static int last_x=0, last_y=0;
static int main_window=0;
static float xy_aspect=1.0;


#define VIS_MAGNETIC_PITCH 100001


class SimulationController
{
public:
  int user_id;
  GLUI_EditText *time_label;
  GLUI_EditText *iter_label;
  GLUI_EditText *imagename_field;
  GLUI_EditText *chkptname_field;
  GLUI_FileBrowser *chkptload_browser;
  enum {
    ACTION_LOAD_FILE,
    ACTION_SAVE_FILE,
    ACTION_TAKE_SCREENSHOT,
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
  void take_screenshot();
  static std::map<int, SimulationController*> instances;
  static void refresh_cb(int user_id);
  static void action_cb(int action_id);
} ;
std::map<int, SimulationController*> SimulationController::instances;


class DatasetController
{
private:
  float DataRange[2];
  GLfloat *ColorData;
  GLfloat *VertexData;
  double AngularOffset;
  int user_id;
  std::vector<int> obj_keys;

  /* live variables */
  int DrawMesh;
  int LogScale;
  int AutoRange;
  int ColorMap;
  int DataMember;
  int Visible;
public:
  DatasetController(GLUI_Panel *parent);
  void SetAngularOffset(double x) { AngularOffset = x; }
  void keyboard(int key, int x, int y) { }
  void draw();
  void refresh();
  double get_scalar(m2vol *V);
  static std::map<int, DatasetController*> instances;
  static void refresh_cb(int user_id);
} ;
std::map<int, DatasetController*> DatasetController::instances;



static SimulationController *sim_controller = NULL;
static m2sim *M2;
static void color_map(double val, GLfloat rgb[3], int ct);
//static void load_next_frame();









void myGlutKeyboard(unsigned char Key, int x, int y)
{
  if (Key == 27) {
    exit(0);
  }
  else {
    sim_controller->keyboard(Key, x, y);
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
  glui->sync_live();
  if (sim_controller->AutoAdvance) {
    m2sim_drive(M2);
    //load_next_frame();
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
  glui = GLUI_Master.create_glui("M2", 0, 768+50, 50);
  new GLUI_StaticText(glui, "M2 controls"); 
  GLUI_Panel *general_panel = new GLUI_Panel(glui, "simulation controls");
  sim_controller = new SimulationController(general_panel);
  GLUI_Panel *dataset_panel = new GLUI_Panel(glui, "data sets");
  DatasetController *ds1 = new DatasetController(dataset_panel);
  new GLUI_Column(dataset_panel);
  DatasetController *ds2 = new DatasetController(dataset_panel);


  ds1->SetAngularOffset(M2_PI);
  ds2->SetAngularOffset(0.0);

  sim_controller->refresh();

  glui->set_main_gfx_window(main_window);
  GLUI_Master.set_glutIdleFunc(myGlutIdle);
  glutMainLoop();
}




SimulationController::SimulationController(GLUI_Panel *parent)
  : rotationX(0.0),
    rotationY(0.0),
    ZoomLevel(-1.0),
    AutoAdvance(0)
{
  user_id = instances.size();
  instances[user_id] = this;
  time_label = new GLUI_EditText(parent, "Time", GLUI_EDITTEXT_FLOAT);
  iter_label = new GLUI_EditText(parent, "Iteration", GLUI_EDITTEXT_INT);
  new GLUI_Checkbox(parent, "Auto advance", &AutoAdvance);
  new GLUI_StaticText(parent, "Zoom level");
  GLUI_Scrollbar *zoom_sb = new GLUI_Scrollbar(parent,
					       "Zoom level",
					       GLUI_SCROLL_HORIZONTAL,
					       &ZoomLevel);
  new GLUI_Separator(parent);
  imagename_field = new GLUI_EditText(parent, "Image file name");
  new GLUI_Button(parent, "Save image", ACTION_TAKE_SCREENSHOT, action_cb);
  new GLUI_Separator(parent);
  new GLUI_Button(parent, "Reset view", ACTION_RESET_VIEW, action_cb);
  new GLUI_Button(parent, "Step", ACTION_STEP, action_cb);
  new GLUI_Button(parent, "Quit", ACTION_QUIT, action_cb);
  new GLUI_Column(parent);
  chkptload_browser = new GLUI_FileBrowser(parent, "Load M2 file");
  new GLUI_Button(parent, "Load checkpoint", ACTION_LOAD_FILE, action_cb);
  chkptname_field = new GLUI_EditText(parent, "Checkpoint file name");
  new GLUI_Button(parent, "Save checkpoint", ACTION_SAVE_FILE, action_cb);
  time_label->deactivate();
  iter_label->deactivate();
  zoom_sb->set_float_limits(-3.0, 0.0);
  imagename_field->set_w(240);
  imagename_field->set_text("m2.ppm");
  chkptname_field->set_w(240);
  chkptname_field->set_text("chkpt.m2");
  chkptload_browser->set_w(200);
  chkptload_browser->set_allow_change_dir(1);
  chkptload_browser->set_show_all_files(0);
  chkptload_browser->add_allowed_extension("m2");
  chkptload_browser->fbreaddir(".");
}
void SimulationController::refresh_cb(int user_id)
{
  instances[user_id]->refresh();
}
void SimulationController::action_cb(int action_id)
{
  switch (action_id) {
  case ACTION_SAVE_FILE:
    m2sim_save_checkpoint(M2, sim_controller->chkptname_field->get_text());
    break;
  case ACTION_LOAD_FILE:
    m2sim_load_checkpoint(M2, sim_controller->chkptload_browser->get_file());
    sim_controller->refresh();
    break;
  case ACTION_TAKE_SCREENSHOT:
    sim_controller->take_screenshot();
    break;
  case ACTION_RESET_VIEW:
    sim_controller->ZoomLevel = -1.0;
    sim_controller->rotationX =  0.0;
    sim_controller->rotationY =  0.0;
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
  time_label->set_float_val(M2->status.time_simulation);
  iter_label->set_int_val(M2->status.iteration_number);
  time_label->deactivate();
  iter_label->deactivate();
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
  glTranslated(0.0, 0.0, -1.0);
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

  char label[256];
  snprintf(label, 256, "t=%3.2f ms", M2->status.time_simulation);
  std::string text = label;
  for (unsigned int i=0; i<text.length(); ++i) {
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, text[i]);
  }
}
void SimulationController::take_screenshot()
{
  int dimx = glutGet(GLUT_WINDOW_WIDTH);
  int dimy = glutGet(GLUT_WINDOW_HEIGHT);
 
  size_t imsize = 3*dimx*dimy;
  char *pixels = (char*) malloc(imsize*sizeof(char));
  glReadPixels(0, 0, dimx, dimy, GL_RGB, GL_UNSIGNED_BYTE, pixels);

  FILE *fp = fopen(imagename_field->get_text(), "wb");
  fprintf(fp, "P6\n%d %d\n255\n", dimx, dimy);
  fwrite(pixels, sizeof(char), imsize, fp);
  fclose(fp);

  free(pixels);
}

DatasetController::DatasetController(GLUI_Panel *parent)
  :
  ColorData(NULL),
  VertexData(NULL),
  AngularOffset(0.0),

  /* live variables */
  DrawMesh(0),
  LogScale(0),
  AutoRange(1),
  ColorMap(3),
  DataMember(0),
  Visible(1)
{
  user_id = instances.size();
  instances[user_id] = this;

  obj_keys.push_back(M2_COMOVING_MASS_DENSITY);
  obj_keys.push_back(M2_GAS_PRESSURE);
  obj_keys.push_back(M2_MAGNETIC_PRESSURE);
  obj_keys.push_back(M2_MACH_NUMBER);
  obj_keys.push_back(M2_MACH_FAST);
  obj_keys.push_back(M2_SIGMA);
  obj_keys.push_back(M2_PLASMA_BETA);
  obj_keys.push_back(M2_ENTROPY);
  obj_keys.push_back(M2_VELOCITY_FOUR_VECTOR0);
  obj_keys.push_back(M2_VELOCITY_FOUR_VECTOR1);
  obj_keys.push_back(M2_VELOCITY_FOUR_VECTOR2);
  obj_keys.push_back(M2_VELOCITY_FOUR_VECTOR3);
  obj_keys.push_back(M2_VELOCITY1);
  obj_keys.push_back(M2_VELOCITY2);
  obj_keys.push_back(M2_VELOCITY3);
  obj_keys.push_back(M2_MAGNETIC1);
  obj_keys.push_back(M2_MAGNETIC2);
  obj_keys.push_back(M2_MAGNETIC3);
  obj_keys.push_back(VIS_MAGNETIC_PITCH);

  GLUI_RadioGroup *radio = new GLUI_RadioGroup(parent, &DataMember,
					       user_id, refresh_cb);
  new GLUI_RadioButton(radio, "mass density");
  new GLUI_RadioButton(radio, "gas pressure");
  new GLUI_RadioButton(radio, "magnetic pressure");
  new GLUI_RadioButton(radio, "sonic mach number");
  new GLUI_RadioButton(radio, "fast mach number");
  new GLUI_RadioButton(radio, "sigma");
  new GLUI_RadioButton(radio, "plasma beta");
  new GLUI_RadioButton(radio, "entropy");
  new GLUI_RadioButton(radio, "u0");
  new GLUI_RadioButton(radio, "u1");
  new GLUI_RadioButton(radio, "u2");
  new GLUI_RadioButton(radio, "u3");
  new GLUI_RadioButton(radio, "v1");
  new GLUI_RadioButton(radio, "v2");
  new GLUI_RadioButton(radio, "v3");
  new GLUI_RadioButton(radio, "B1");
  new GLUI_RadioButton(radio, "B2");
  new GLUI_RadioButton(radio, "B3");
  new GLUI_RadioButton(radio, "B-pitch");
  new GLUI_Separator(parent);
  new GLUI_Checkbox(parent, "draw mesh", &DrawMesh, user_id, refresh_cb);
  new GLUI_Checkbox(parent, "log scale", &LogScale, user_id, refresh_cb);
  new GLUI_Checkbox(parent, "auto range", &AutoRange, user_id, refresh_cb);
  new GLUI_Checkbox(parent, "visible", &Visible, user_id, refresh_cb);
  GLUI_Spinner *cm = new GLUI_Spinner (parent, "color map", &ColorMap, user_id, refresh_cb);
  new GLUI_Spinner (parent, "range min",
		    &DataRange[0], user_id, refresh_cb);
  new GLUI_Spinner (parent, "range max",
		    &DataRange[1], user_id, refresh_cb);
  cm->set_int_limits(0, 5);
  this->refresh();
}
void DatasetController::draw()
{
  if (!Visible) return;
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
double DatasetController::get_scalar(m2vol *V)
{
  int mem = obj_keys[DataMember];
  double rhat[4] = { 0.0, 1.0, 0.0, 0.0 };
  double y;
  if (mem == VIS_MAGNETIC_PITCH) {
    y = atan(m2aux_get(&V->aux, M2_MAGNETIC3)/
	     m2aux_get(&V->aux, M2_MAGNETIC1)) * 180 / M2_PI;
  }
  else if (mem == M2_MACH_ALFVEN ||
      mem == M2_MACH_FAST ||
      mem == M2_MACH_SLOW) {
    y = m2aux_mach(&V->aux, rhat, mem);
  }
  else {
    y = m2aux_get(&V->aux, mem);
  }
  if (LogScale) {
    if (fabs(y) < 1e-6) {
      y = 1e-6;
    }
    y = log10(fabs(y));
  }
  return y;
}
void DatasetController::refresh()
{
  if (!Visible) return;
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
  if (AutoRange) {
    DataRange[0] = +1e10;
    DataRange[1] = -1e10;
    for (int n=0; n<L[0]; ++n) {
      V = M2->volumes + n;
      if (V->zone_type != M2_ZONE_TYPE_FULL) {
	continue;
      }
      y = get_scalar(V);
      if (y < DataRange[0]) DataRange[0] = y;
      if (y > DataRange[1]) DataRange[1] = y;
    }
  }
  for (int n=0; n<L[0]; ++n) {
    V = M2->volumes + n;
    if (V->zone_type != M2_ZONE_TYPE_FULL) {
      VertexData[4*3*n + 0*3 + 0] = 0.0;
      VertexData[4*3*n + 0*3 + 1] = 0.0;
      VertexData[4*3*n + 0*3 + 2] = 0.0;
      VertexData[4*3*n + 1*3 + 0] = 0.0;
      VertexData[4*3*n + 1*3 + 1] = 0.0;
      VertexData[4*3*n + 1*3 + 2] = 0.0;
      VertexData[4*3*n + 2*3 + 0] = 0.0;
      VertexData[4*3*n + 2*3 + 1] = 0.0;
      VertexData[4*3*n + 2*3 + 2] = 0.0;
      VertexData[4*3*n + 3*3 + 0] = 0.0;
      VertexData[4*3*n + 3*3 + 1] = 0.0;
      VertexData[4*3*n + 3*3 + 2] = 0.0;
      ColorData[3*n + 0] = 0.0;
      ColorData[3*n + 1] = 0.0;
      ColorData[3*n + 2] = 0.0;
      continue;
    }
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
    y = get_scalar(M2->volumes + n);
    y -= DataRange[0];
    y /= DataRange[1] - DataRange[0];
    color_map(y, &ColorData[3*n], ColorMap);
  }
}
void DatasetController::refresh_cb(int user_id)
{
  instances[user_id]->refresh();
}



void load_next_frame()
{
  static int checkpoint_number = 0;
  char fname[256];
  if (checkpoint_number == 900) {
    return;
  }
  else {
    checkpoint_number += 1;
  }
  snprintf(fname, 256, "sgrbA/chkpt.%04d.m2", checkpoint_number);
  m2sim_load_checkpoint(M2, fname);
  sim_controller->refresh();

  if (M2->status.time_simulation < 100.0) {
    sim_controller->ZoomLevel -= 0.001;
  }
  else if (M2->status.time_simulation < 1000.0) {
    if (sim_controller->ZoomLevel > -2.5) {
      sim_controller->ZoomLevel -= 0.005;
    }
  }

  snprintf(fname, 256, "sgrbA/pressure.%04d.ppm", checkpoint_number);
  sim_controller->imagename_field->set_text(fname);
  sim_controller->take_screenshot();
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
