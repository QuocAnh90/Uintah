#include <Core/Datatypes/GLVolRenState.h>
#include <Core/Datatypes/GLVolumeRenderer.h>
#include <Core/Datatypes/GLTexture3D.h>
#include <Core/Datatypes/Brick.h>
#include <Core/Datatypes/Polygon.h>
#include <GL/gl.h>
#include <vector>
#include <stdlib.h>
#include <iostream>

namespace SCIRun {

using std::vector;
using std::cerr;
using std::endl;


GLVolRenState::GLVolRenState(const GLVolumeRenderer* glvr)
    : volren( glvr ), texName(0), reload((unsigned char *)1)
{
  // Base Class, holds pointer to VolumeRenderer and 
  // common computation
}


void
GLVolRenState::computeView(Ray& ray)
{
  double mvmat[16];
  Transform mat;
  Vector view;
  Point viewPt;
      
  glGetDoublev( GL_MODELVIEW_MATRIX, mvmat);
  /* remember that the glmatrix is stored as
       0  4  8 12
       1  5  9 13
       2  6 10 14
       3  7 11 15 */
 
  view = Vector(mvmat[12], mvmat[13], mvmat[14]);
  view.normalize();
  viewPt = Point(-mvmat[12], -mvmat[13], -mvmat[14]);
    
  /* set the translation to zero */
  mvmat[12] = mvmat[13] = mvmat[14] = 0;
  /* Because of the order of the glmatrix we are storing as a transpose.
       if there is not use of scale then the transpose is the  inverse */
  mat.set( mvmat );
    
  /* project view info into object space */
  view = mat.project( view );
  viewPt = mat.project( viewPt );


  ray =  Ray(viewPt, view);
}

void
GLVolRenState::drawPolys( vector<Polygon *> polys )
{
  int i;
  
  for (i = 0; i < polys.size(); i++) {
    switch (polys[i]->size() ) {
    case 1:
      glBegin(GL_POINTS);
      glVertex3f((*(polys[i]))[0].x(),(*(polys[i]))[0].y(),
		 (*(polys[i]))[0].z());
      glEnd();
      break;
    case 2:
      glBegin(GL_LINES);
      glVertex3f((*(polys[i]))[0].x(),(*(polys[i]))[0].y(),
		 (*(polys[i]))[0].z());
      glVertex3f((*(polys[i]))[1].x(), (*(polys[i]))[1].y(),
		 (*(polys[i]))[1].z());
      glEnd();
      break;
    case 3:
      glBegin(GL_TRIANGLES);
      glVertex3f((*(polys[i]))[0].x(),(*(polys[i]))[0].y(),
		 (*(polys[i]))[0].z());
      glVertex3f((*(polys[i]))[1].x(), (*(polys[i]))[1].y(),
		 (*(polys[i]))[1].z());
      glVertex3f((*(polys[i]))[2].x(),(*(polys[i]))[2].y(),
		 (*(polys[i]))[2].z());
      glEnd();
      break;
    case 4:
    case 5:
    case 6:
      {
	int k;
	glBegin(GL_POLYGON);
	for(k =0; k < polys[i]->size(); k++)
	{
	  glVertex3f((*(polys[i]))[k].x(),(*(polys[i]))[k].y(),
		     (*(polys[i]))[k].z());
	}
	glEnd();
      }
      break;
    }
  }
}

void
GLVolRenState::loadColorMap(Brick& brick)
{
#ifdef __sgi
  glColorTable(GL_TEXTURE_COLOR_TABLE_SGI,
               GL_RGBA,
               256, // try larger sizes?
               GL_RGBA,  // need an alpha value...
               GL_UNSIGNED_BYTE, // try shorts...
               volren->TransferFunctions[brick.level()]);
#endif
}

void 
GLVolRenState::loadTexture(Brick& brick)
{
  if( !brick.texName() || reload ) {
    if( !brick.texName() )
      glGenTextures(1, brick.texNameP());

    glBindTexture(GL_TEXTURE_3D_EXT, brick.texName());

    if( volren->_interp ){
      glTexParameteri(GL_TEXTURE_3D_EXT, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_3D_EXT, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    } else {
      glTexParameteri(GL_TEXTURE_3D_EXT, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_3D_EXT, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    }

    glTexParameteri(GL_TEXTURE_3D_EXT, GL_TEXTURE_WRAP_S,
		    GL_CLAMP);
    glTexParameteri(GL_TEXTURE_3D_EXT, GL_TEXTURE_WRAP_T,
		    GL_CLAMP);
    glTexParameteri(GL_TEXTURE_3D_EXT, GL_TEXTURE_WRAP_R_EXT,
		    GL_CLAMP);

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    
    // set up the texture
    //glTexImage3DEXT(GL_TEXTURE_3D_EXT, 0,
    glTexImage3D(GL_TEXTURE_3D_EXT, 0,
		    GL_INTENSITY8,
		    (brick.texture())->dim1(), 
		    (brick.texture())->dim2(), 
		    (brick.texture())->dim3(),
		    0,
		    GL_RED, GL_UNSIGNED_BYTE,
		    &(*(brick.texture()))(0,0,0));

  } else {
    glBindTexture(GL_TEXTURE_3D_EXT, brick.texName());
  }
}
void 
GLVolRenState::makeTextureMatrix( const Brick& brick)
{
  double splane[4]={0,0,0,0};
  double tplane[4]={0,0,0,0};
  double rplane[4]={0,0,0,0};
  double qplane[4]={0,0,0,1};


  Vector diag;

  
  /* The cube is numbered in the following way 
      
         2________ 6        y
        /|       /|         |  
       / |      / |         |
      /  |     /  |         |
    3/__0|____/7__|4        |_________ x
     |   /    |   /         /
     |  /     |  /         /
     | /      | /         /
    1|/_______|/5        /
                        z  
  */



  diag = brick[7] - brick[0];


  glTexGend(GL_S,GL_TEXTURE_GEN_MODE,GL_OBJECT_LINEAR);
  glTexGend(GL_T,GL_TEXTURE_GEN_MODE,GL_OBJECT_LINEAR);
  glTexGend(GL_R,GL_TEXTURE_GEN_MODE,GL_OBJECT_LINEAR);
  glTexGend(GL_Q,GL_TEXTURE_GEN_MODE,GL_OBJECT_LINEAR);

  //  This code is for render overlapping bricks.  The plane equations
  //  for s are  (Nx * Pxmin) + d = aX/2  and
  //  (Nx * Pxmax) + d = 1 - aX/2 where
  //  Nx is the x component of the normal,  Pxmin and Pxmax are the x 
  //  components of the min and max points on the TexCube, and  aX is one
  //  texel width.  Solving for Nx and d we get
  //  Nx = (1 - aX)/(Pxmax - Pxmin) and
  //  d = aX/2 - (Pxmin *(1 - aX))/(Pxmax - Pxmin)

  splane[0] = (1 - brick.aX * (brick.padx + 1))/diag.x();
  splane[3] = brick.aX * 0.5 - (brick[0].x() *
				(1 - brick.aX * (brick.padx+1))/diag.x());
  tplane[1] = (1 - brick.aY * (brick.pady + 1))/diag.y();
  tplane[3] = brick.aY * 0.5 - (brick[0].y() *
				(1 - brick.aY * (brick.pady+1))/diag.y());
  rplane[2] = (1 - brick.aZ * (brick.padz + 1))/diag.z();
  rplane[3] = brick.aZ * 0.5 - (brick[0].z() *
				(1 - brick.aZ * (brick.padz+1))/diag.z());

  
  glTexGendv(GL_S,GL_OBJECT_PLANE,splane);
  glTexGendv(GL_T,GL_OBJECT_PLANE,tplane);
  glTexGendv(GL_R,GL_OBJECT_PLANE,rplane);
  glTexGendv(GL_Q,GL_OBJECT_PLANE,qplane);
}

void 
GLVolRenState::enableTexCoords()
{
  glEnable(GL_TEXTURE_GEN_S);
  glEnable(GL_TEXTURE_GEN_T);
  glEnable(GL_TEXTURE_GEN_R);
  glEnable(GL_TEXTURE_GEN_Q);
}
void 
GLVolRenState::disableTexCoords()
{
  glDisable(GL_TEXTURE_GEN_S);
  glDisable(GL_TEXTURE_GEN_T);
  glDisable(GL_TEXTURE_GEN_R);
  glDisable(GL_TEXTURE_GEN_Q);
}

void 
GLVolRenState::enableBlend()
{
  glEnable(GL_BLEND);
}
void 
GLVolRenState::disableBlend()
{
  glDisable(GL_BLEND);
}

void
GLVolRenState::drawWireFrame(const Brick& brick)
{
  int i;
  glEnable(GL_DEPTH_TEST);
//   double r,g,b;
//   char c;
//   r = drand48();
//   g = drand48();
//   b = drand48();
//   std::cin.get(c);
//   glColor4f(r,g,b,1.0);
  glColor4f(0.8,0.8,0.8,1.0);
  glPushMatrix();
  glBegin(GL_LINES);
  for(i = 0; i < 4; i++){
    glVertex3d(brick[i].x(), brick[i].y(), brick[i].z());
    glVertex3d(brick[i+4].x(), brick[i+4].y(), brick[i+4].z());
  }
  glEnd();

  glBegin(GL_LINE_LOOP);
   glVertex3d(brick[0].x(), brick[0].y(), brick[0].z());
   glVertex3d(brick[1].x(), brick[1].y(), brick[1].z());
   glVertex3d(brick[3].x(), brick[3].y(), brick[3].z());
   glVertex3d(brick[2].x(), brick[2].y(), brick[2].z());
  glEnd();

  glBegin(GL_LINE_LOOP);
   glVertex3d(brick[4].x(), brick[4].y(), brick[4].z());
   glVertex3d(brick[5].x(), brick[5].y(), brick[5].z());
   glVertex3d(brick[7].x(), brick[7].y(), brick[7].z());
   glVertex3d(brick[6].x(), brick[6].y(), brick[6].z());
  glEnd();
  glPopMatrix();
  glDisable(GL_DEPTH_TEST);

}

} // End namespace SCIRun


