//  
//  For more information, please see: http://software.sci.utah.edu
//  
//  The MIT License
//  
//  Copyright (c) 2004 Scientific Computing and Imaging Institute,
//  University of Utah.
//  
//  
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//  
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//  
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//  
//    File   : TextureRenderer.h
//    Author : Milan Ikits
//    Date   : Wed Jul  7 23:34:33 2004

#ifndef TextureRenderer_h
#define TextureRenderer_h

#include <Core/Thread/Mutex.h>
#include <Core/Geom/ColorMap.h>
#include <Core/Geom/GeomObj.h>

#include <Core/Containers/Array2.h>
#include <Core/Containers/Array3.h>

#include <Core/Volume/TextureBrick.h>
#include <Core/Volume/Texture.h>
#include <Core/Volume/ColorMap2.h>
#include <Core/Volume/CM2Widget.h>

#include <Core/Geometry/Polygon.h>

#include <Core/Volume/share.h>

namespace SCIRun {

class Pbuffer;
class FragmentProgramARB;
class VolShaderFactory;

class SCISHARE TextureRenderer : public GeomObj
{
public:
  TextureRenderer(TextureHandle tex, 
		  ColorMapHandle cmap1, 
		  vector<ColorMap2Handle> & cmap2,
                  int tex_mem);
  TextureRenderer(const TextureRenderer&);
  virtual ~TextureRenderer();

  void set_texture(TextureHandle tex);
  void set_colormap1(ColorMapHandle cmap1);
  void set_colormap2(vector<ColorMap2Handle> &cmap2);
  void set_colormap_size(int size);
  void set_slice_alpha(double alpha);
  void set_sw_raster(bool b);
  bool use_pbuffer() { return use_pbuffer_; }
  void set_blend_num_bits(int b);
  bool use_blend_buffer();
  void set_stencil(bool use){ use_stencil_ = use; }
  void invert_opacity(bool invert){ invert_opacity_ = invert; }
  inline void set_interp(bool i) { interp_ = i; }
  
  enum RenderMode { MODE_NONE, MODE_OVER, MODE_MIP, MODE_SLICE };
    
#ifdef SCI_OPENGL
  virtual void draw(DrawInfoOpenGL*, Material*, double time) = 0;
#endif
  
  bool saveobj(std::ostream&, const std::string& format, GeomSave*);
  virtual GeomObj* clone() = 0;
  virtual void get_bounds(BBox& bb) { tex_->get_bounds(bb); }
  virtual void io(Piostream&);
  static PersistentTypeID type_id;

protected:
  TextureHandle tex_;
  Mutex mutex_;
  ColorMapHandle cmap1_;
  vector<ColorMap2Handle> &cmap2_;
  bool cmap1_dirty_;
  bool cmap2_dirty_;
  bool alpha_dirty_;
  RenderMode mode_;
  bool interp_;
  int lighting_;
  double sampling_rate_;
  double irate_;
  bool imode_;
  double slice_alpha_;
  bool hardware_raster_;
  DrawInfoOpenGL* di_;

#ifdef SCI_OPENGL  
  int cmap2_size_;
  Array2<float> cmap1_array_;
  unsigned int cmap1_tex_;
  Array3<float> raster_array_;
  Array3<unsigned char> cmap2_array_;
  unsigned int cmap2_tex_;
  bool use_pbuffer_;
  Pbuffer* raster_pbuffer_;
  CM2ShaderFactory* shader_factory_;
  Pbuffer* cmap2_pbuffer_;
  FragmentProgramARB* cmap2_shader_nv_;
  FragmentProgramARB* cmap2_shader_ati_;
  VolShaderFactory* vol_shader_factory_;
  Pbuffer* blend_buffer_;
  int blend_num_bits_;
  bool use_blend_buffer_;
  int free_tex_mem_;
  bool use_stencil_;
  bool invert_opacity_;
  
  struct TexParam
  {
    int nx, ny, nz, nb;
    uint id;
    TextureBrickHandle brick;
    int comp;
    GLenum textype;
    TexParam() : nx(0), ny(0), nz(0), nb(0), id(0), brick(0), comp(0), textype(GL_UNSIGNED_BYTE) {}
    TexParam(int x, int y, int z, int b, GLenum f, uint i)
      : nx(x), ny(y), nz(z), nb(b), id(i), brick(0), comp(0), textype(f) {}
  };
  vector<TexParam> tex_pool_;
  
  // Tests the bounding box against the current MODELVIEW and
  // PROJECTION matrices to determine if it is within the viewport.
  // Returns true if it is visible.
  bool test_against_view(const BBox &bbox);

  Ray compute_view();
  void load_brick(vector<TextureBrickHandle> &b, int i, bool use_cmap2);
  void draw_polygons(vector<float>& vertex, vector<float>& texcoord,
		     vector<int>& poly,
                     bool normal, bool fog, Pbuffer* buffer,
		     vector<int> *mask = 0, FragmentProgramARB *shader = 0);
  void draw_polygons_wireframe(vector<float>& vertex, vector<float>& texcoord,
			       vector<int>& poly,
			       bool normal, bool fog, Pbuffer* buffer,
			       vector<int> *mask=0);

  void build_colormap1(Array2<float>& cmap_array,
		       unsigned int& cmap_tex, bool& cmap_dirty,
		       bool& alpha_dirty,  double level_exponent = 0.0);

  void build_colormap2();
  void colormap2_hardware_rasterize_setup();
  void colormap2_hardware_destroy_buffers();
  void colormap2_hardware_rasterize();
  void colormap2_software_rasterize();



  void bind_colormap1(Array2<float>& cmap_array, unsigned int cmap_tex );
  void bind_colormap2();

  void release_colormap1();
  void release_colormap2();

#endif
};



} // end namespace SCIRun

#endif // TextureRenderer_h
