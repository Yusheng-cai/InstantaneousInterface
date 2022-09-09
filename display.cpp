#include "src/Mesh.h"
#include "tools/CommandLineArguments.h"
#include "tools/CommonTypes.h"

#include <GL/glut.h>
# include <GL/freeglut.h>
# include <GL/gl.h>
#include <vector>
#include <array>

#  define GLUT_WHEEL_UP   3
#  define GLUT_WHEEL_DOWN 4

using Real3 = CommonTypes::Real3;
using INT3  = CommonTypes::index3;

Mesh m;
float _table_angle = 0.0f;
bool _3d_view = true;
static int _win_number;


void draw_mesh(Mesh& mesh);
void display();
void key_stroke (unsigned char c, int mouseX, int mouseY);
void mouse_keys (int button, int state, int x, int y);


int main(int argc, char** argv)
{
    std::string fname = argv[1];
    MeshTools::readPLYlibr(fname, m);

    glutInit (&argc, argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize (900, 900);
    glutInitWindowPosition (240, 212);
    _win_number = glutCreateWindow (argv[0]);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glLineWidth(2.2f);
    //glEnable(GL_LINE_SMOOTH);
    //glLightfv(GL_LIGHT0, GL_AMBIENT,  {});
    //glLightfv(GL_LIGHT0, GL_DIFFUSE,  Global::_light0_diffuse);
    //GLfloat lightp[] = {1.0f,1.0f,1.0f,1.0f};
    //glLightfv(GL_LIGHT0, GL_SPECULAR, lightp);

    glutKeyboardFunc(key_stroke);
    glutMouseFunc(mouse_keys);
    glutDisplayFunc(display);
    glutMainLoop();

    return 0;
}


void draw_mesh(Mesh& mesh) 
{
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);

    mesh.CalcVertexNormals();
    const auto& v = mesh.getvertices();
    const auto& t = mesh.gettriangles();
    std::vector<Real3> verts(v.size());
    std::vector<INT3>  tris(t.size());
    std::vector<Real3> norms(v.size());

    for (int i=0;i<v.size();i++)
    {
        verts[i] = v[i].position_;
        norms[i] = v[i].normals_;
    }

    for (int i=0;i<t.size();i++)
    {
        tris[i] = t[i].triangleindices_;
    }

    glVertexPointer(3, GL_FLOAT, 0, verts.data());
    glNormalPointer(GL_FLOAT, 0, norms.data());
    //glColorPointer(3, GL_FLOAT, 0, mesh._colors.data());
    glDrawElements(GL_TRIANGLES, tris.size()*3, GL_UNSIGNED_INT, tris.data());

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
}

void display()
{
    glClearColor ( 0.93,  0.93,  0.93,  0.93);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glViewport(0,0,900,900);

    glMatrixMode (GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60., 1., 0.5, 100.);

    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity ();

    if(_3d_view){
        GLfloat lightPos0[] = {0.0f, 0.0f, 0.0f, 1.0f};
        glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);

        glTranslatef(0.0, 0.0, -1.5);
        glRotatef(45.0f, -1.0f, 0.0f, 0.0f);
        glRotatef(_table_angle, -1.0f, 0.0f, 0.0f);

        static float angle = 0.0f;
        angle = fmodf(angle+0.3f, 360.f);
        //glRotatef(angle, 1.0f, 0.0f, 0.0f);
        glRotatef(angle, 0.0f, 0.0f, 1.0f);
        glutPostRedisplay();
    }
    else
        glTranslatef(0.0, 0.0, -1.0);
    float s = 0.5f;
    glScalef(s, s, s);

    draw_mesh(m);

    glutSwapBuffers();
    glFlush ();
}

void key_stroke (unsigned char c, int mouseX, int mouseY) {
    static bool wires  = false;

    switch (c) {
    case 27 :
        glFinish();
        glutDestroyWindow(_win_number);
        exit (0);
        break;
    case 'w' :
        if(!wires) {
            glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
            glDisable(GL_LIGHTING);
            glutPostRedisplay();
            wires = true;
        }else {
            glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
            glEnable(GL_LIGHTING);
            glutPostRedisplay();
            wires = false;
        }
        break;
    }
}


void mouse_keys (int button, int state, int x, int y)
{
    if(button == GLUT_WHEEL_UP)
        _table_angle += 1.0f;
    if(button == GLUT_WHEEL_DOWN)
        _table_angle -= 1.0f;
}