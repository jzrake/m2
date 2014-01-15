extern "C" {
#include "m2.h" 
}
#if (M2_HAVE_GLUT)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <GL/glut.h>
#define ESCAPE_KEY 27
 

/* ------------------------------------------ */
/* GLUT callback functions and object handles */
/* ------------------------------------------ */
static void GLUTDisplayFunc();
static void GLUTIdleFunc();
static void GLUTReshapeFunc(int width, int height);
static void GLUTKeyboardFunc(unsigned char key, int x, int y);
static void GLUTSpecialFunc(int key, int x, int y);

static int GlutWindow;



/* ----------------------------------- */
/* internal utility functions and data */
/* ----------------------------------- */
static void menu_callback_main(int num);
static void menu_callback_field_variable(int num);
static void menu_callback_color_mapping(int num);
static void create_menu();
static void reload_rgb_data();
static void color_map(double val, GLfloat rgb[3]);

static double ZoomLevel = 1.0;
static double TranslateZ = -2.0;
static double RotationX = 0.0;
static double RotationY = 0.0;
static m2sim *M2; 

static int DataMember = M2_COMOVING_MASS_DENSITY;
static int ColorTable = 3;
static int LogScaling = 0;
static int AutoPlay = 0;
static int DrawMesh = 0;
static double DataRange[2];
static GLfloat *VertexData = NULL;
static GLfloat *ColorData = NULL;

enum {
  VIS_EXIT,
  VIS_TOGGLE_LOG_SCALE,
  VIS_TOGGLE_DRAW_MESH,
  VIS_TOGGLE_AUTO_PLAY,
} ;


void m2sim_visualize(m2sim *m2, int argc, char **argv)
{
  M2 = m2;
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);
  glutInitWindowSize(768, 768);
  glutInitWindowPosition(512, 0);
  GlutWindow = glutCreateWindow("m2");
  glutDisplayFunc(GLUTDisplayFunc);
  glutIdleFunc(GLUTIdleFunc);
  glutReshapeFunc(GLUTReshapeFunc);
  glutKeyboardFunc(GLUTKeyboardFunc);
  glutSpecialFunc(GLUTSpecialFunc);

  create_menu();
  reload_rgb_data();

  glutMainLoop();
}


/* --------------------------------------------------------------------------
 * GLUT function callbacks
 * --------------------------------------------------------------------------
 */
void GLUTDisplayFunc()
{
  glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
 
  /* reset and configure camera */
  glLoadIdentity();
  glTranslated(0.0, 0.0, TranslateZ);
  //  glRotated(RotationZ, 0.0, 0.0, 1.0);
  glRotated(RotationX, 1.0, 0.0, 0.0);
  glRotated(RotationY, 0.0, 1.0, 0.0);
  glRotated(-90.0, 1.0, 0.0, 0.0);
  glScaled(ZoomLevel, ZoomLevel, ZoomLevel);
 
  /* actually draw stuff */
  glEnable(GL_DEPTH_TEST);


  int n;
  m2vol *V;
  for (n=0; n<M2->local_grid_size[0]; ++n) {
    V = M2->volumes + n;

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

 
  /* end drawing */
  glutSwapBuffers();

  if (AutoPlay) {
    m2sim_drive(M2);
    reload_rgb_data();
  }
}
 
void GLUTIdleFunc()
{
  glutPostRedisplay();
}
 
void GLUTReshapeFunc(int width, int height)
{
  if (height == 0) {
    height = 1;
  }
  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0f, (GLfloat)width/(GLfloat)height, 0.1f, 100.0f);
  glMatrixMode(GL_MODELVIEW);
}
 
void GLUTKeyboardFunc(unsigned char key, int x, int y)
{
  switch (key) {
  case ESCAPE_KEY: exit(0);
  case '-': ZoomLevel /= 1.1; break;
  case '=': ZoomLevel *= 1.1; break;
  case 's': m2sim_drive(M2); reload_rgb_data(); break;
  case 'p': printf("data range: [%4.3e %4.3e]\n", DataRange[0], DataRange[1]);
    break;
  }
}
 
void GLUTSpecialFunc(int key, int x, int y)
{
  double a = 4.0;
  switch (key) {
  case GLUT_KEY_RIGHT: RotationY -= a; break;
  case GLUT_KEY_LEFT:  RotationY += a; break;
  case GLUT_KEY_UP:    RotationX -= a; break;
  case GLUT_KEY_DOWN:  RotationX += a; break;
  }
  glutPostRedisplay();
}




/* private module functions */
void menu_callback_field_variable(int num)
{
  DataMember = num;
  reload_rgb_data();
  glutPostRedisplay();
}
void menu_callback_color_mapping(int num)
{
  ColorTable = num;
  reload_rgb_data();
  glutPostRedisplay();
}
void menu_callback_main(int num)
{
  switch (num) {
  case VIS_EXIT:
    glutDestroyWindow(GlutWindow);
    exit(0);
    break;
  case VIS_TOGGLE_LOG_SCALE:
    LogScaling ^= 1;
    reload_rgb_data();
    printf("Log scaling is %s\n", LogScaling ? "enabled" : "disabled");
    break;
  case VIS_TOGGLE_DRAW_MESH:
    DrawMesh ^= 1;
    break;
  case VIS_TOGGLE_AUTO_PLAY:
    AutoPlay ^= 1;
    break;
  }
}
void create_menu()
{
  int main_menu;
  int field_variable;
  int color_mapping;


  field_variable = glutCreateMenu(menu_callback_field_variable);
  glutAddMenuEntry("velocity four vector: 0", M2_VELOCITY_FOUR_VECTOR0);
  glutAddMenuEntry("velocity four vector: 1", M2_VELOCITY_FOUR_VECTOR1);
  glutAddMenuEntry("velocity four vector: 2", M2_VELOCITY_FOUR_VECTOR2);
  glutAddMenuEntry("velocity four vector: 3", M2_VELOCITY_FOUR_VECTOR3);
  glutAddMenuEntry("magnetic four vector: 0", M2_MAGNETIC_FOUR_VECTOR0);
  glutAddMenuEntry("magnetic four vector: 1", M2_MAGNETIC_FOUR_VECTOR1);
  glutAddMenuEntry("magnetic four vector: 2", M2_MAGNETIC_FOUR_VECTOR2);
  glutAddMenuEntry("magnetic four vector: 3", M2_MAGNETIC_FOUR_VECTOR3);
  glutAddMenuEntry("comoving mass density", M2_COMOVING_MASS_DENSITY);
  glutAddMenuEntry("gas pressure", M2_GAS_PRESSURE);
  glutAddMenuEntry("magnetic pressure", M2_MAGNETIC_PRESSURE);
  glutAddMenuEntry("sigma", M2_SIGMA);


  color_mapping = glutCreateMenu(menu_callback_color_mapping);
  glutAddMenuEntry("color table 0", 0);
  glutAddMenuEntry("color table 1", 1);
  glutAddMenuEntry("color table 2", 2);
  glutAddMenuEntry("color table 3", 3);
  glutAddMenuEntry("color table 4", 4);
  glutAddMenuEntry("color table 5", 5);


  main_menu = glutCreateMenu(menu_callback_main);
  glutAddSubMenu("Field Variable", field_variable);
  glutAddSubMenu("Color Mapping", color_mapping);
  glutAddMenuEntry("Toggle log scale", VIS_TOGGLE_LOG_SCALE);
  glutAddMenuEntry("Toggle draw mesh", VIS_TOGGLE_DRAW_MESH);
  glutAddMenuEntry("Toggle auto play", VIS_TOGGLE_AUTO_PLAY);
  glutAddMenuEntry("Quit", VIS_EXIT);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
}
void reload_rgb_data()
{
  int n;
  int *L = M2->local_grid_size;
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

  VertexData = (GLfloat*) realloc(VertexData, L[0] * 3 * 4 * sizeof(GLfloat));
  ColorData = (GLfloat*) realloc(ColorData, L[0] * 3 * sizeof(GLfloat));

  for (n=0; n<L[0]; ++n) {
    V = M2->volumes + n;

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
    V = M2->volumes + n;

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

    // x00[1] = V->x0[1];
    // x00[2] = V->x0[2];
    // x00[3] = V->x0[3];
    // x01[1] = V->x1[1];
    // x01[2] = V->x0[2];
    // x01[3] = V->x0[3];
    // x11[1] = V->x1[1];
    // x11[2] = V->x0[2];
    // x11[3] = V->x1[3];
    // x10[1] = V->x0[1];
    // x10[2] = V->x0[2];
    // x10[3] = V->x1[3];

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

    y = m2aux_get(&V->aux, DataMember);

    if (LogScaling) {
      y = log10(fabs(y));
    }

    y -= DataRange[0];
    y /= DataRange[1] - DataRange[0];
    color_map(y, &ColorData[3*n]);
  }
}
void color_map(double val, GLfloat rgb[3])
{
  double rrr, ggg, bbb;
  if( ColorTable == 0 ){
    double nexp = 8.0;
    rrr = exp(-nexp*pow(val-5./6.,2.0)) + .25*exp(-nexp*pow(val+1./6.,2.0));
    ggg = exp(-nexp*pow(val-3./6.,2.0));
    bbb = exp(-nexp*pow(val-1./6.,2.0)) + .25*exp(-nexp*pow(val-7./6.,2.0));
  }else if(ColorTable == 1){
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
  }else if(ColorTable == 2){
    rrr = 2.*val;
    ggg = 1.2*val;
    bbb = .8*val;
  }else if(ColorTable == 3){
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
  }else if(ColorTable == 4){
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
  }else if(ColorTable == 5){
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
