#include <stdlib.h>
#include <stdio.h>
#include <GL/glut.h>
#define ESCAPE_KEY 27
#include "m2.h" 
 
static void GLUTDisplayFunc();
static void GLUTIdleFunc();
static void GLUTReshapeFunc(int width, int height);
static void GLUTKeyboardFunc(unsigned char key, int x, int y);
static void GLUTSpecialFunc(int key, int x, int y);
 
static int GlutWindow;
static double TranslateZ = -2.0;
static double RotationX = 90.0;
static double RotationY =  0.0;
//static double RotationZ = 0.0;
static m2sim *M2; 

 
static int menu_id;
static int submenu_id;
static int value = 0;
static void menu_callback(int num)
{
  if (num == 0) {
    glutDestroyWindow(GlutWindow);
    exit(0);
  }
  else {
    value = num;
  }
  glutPostRedisplay();
}
static void create_menu()
{
  submenu_id = glutCreateMenu(menu_callback);
  glutAddMenuEntry("Sphere", 2);
  glutAddMenuEntry("Cone", 3);
  glutAddMenuEntry("Torus", 4);
  glutAddMenuEntry("Teapot", 5);
  menu_id = glutCreateMenu(menu_callback);
  glutAddSubMenu("Draw", submenu_id);
  glutAddMenuEntry("Clear", 1);
  glutAddMenuEntry("Quit", 0);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
}


void m2sim_visualize(m2sim *m2, int argc, char **argv)
{
  M2 = m2;
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);
  glutInitWindowSize(768, 768);
  glutInitWindowPosition(512, 0);
  GlutWindow = glutCreateWindow("m2");
  create_menu();
  glutDisplayFunc(GLUTDisplayFunc);
  glutIdleFunc(GLUTIdleFunc);
  glutReshapeFunc(GLUTReshapeFunc);
  glutKeyboardFunc(GLUTKeyboardFunc);
  glutSpecialFunc(GLUTSpecialFunc);
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
  glTranslated(0.0, 0.0, 0.0);
  glTranslated(0.0, 0.0, TranslateZ);
  //  glRotated(RotationZ, 0.0, 0.0, 1.0);
  glRotated(RotationX, 1.0, 0.0, 0.0);
  glRotated(RotationY, 0.0, 1.0, 0.0);
  glRotated(-90.0, 1.0, 0.0, 0.0);
 
  /* actually draw stuff */
  glEnable(GL_DEPTH_TEST);


  int n;
  m2vol *V;
  for (n=0; n<M2->local_grid_size[0]; ++n) {
    V = M2->volumes + n;
    //glBegin(GL_LINE_LOOP);
    glBegin(GL_QUADS);
    glColor3d(V->prim.B3, V->prim.B3, V->prim.B3);
    glVertex3d(V->x0[1], V->x0[2], 0.0);
    glVertex3d(V->x1[1], V->x0[2], 0.0);
    glVertex3d(V->x1[1], V->x1[2], 0.0);
    glVertex3d(V->x0[1], V->x1[2], 0.0);
    glEnd();
  }

 
  /* end drawing */
  glutSwapBuffers();
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
  case '-': TranslateZ *= 1.1; break;
  case '=': TranslateZ /= 1.1; break;
  case 's': m2sim_drive(M2); break;
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

