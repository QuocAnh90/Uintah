/*
 *  Salmon.cc:  The Geometry Viewer Window
 *
 *  Written by:
 *   David Weinstein
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 *  Copyright (C) 1994 SCI Group
 */

// Someday, we should delete these four lines, when the
// compiler stops griping about const cast away...
#include <X11/Intrinsic.h>
#include "myStringDefs.h"
#include "myXmStrDefs.h"
#include "myShell.h"

#include <Geom.h>
#include <XQColor.h>
#include <Salmon/Salmon.h>
#include <Salmon/Roe.h>
#include <MotifCallback.h>
#include <MtXEventLoop.h>
#include <NetworkEditor.h>
#include <NotFinished.h>
#include <Mt/DialogShell.h>
#include <Mt/DrawingArea.h>
#include <Mt/Form.h>
#include <Mt/Frame.h>
#include <Mt/GLwMDraw.h>
#include <Mt/RowColumn.h>
#include <Mt/Label.h>
#include <Mt/ScrolledWindow.h>
#include <Mt/ToggleButton.h>
#include <Mt/PushButton.h>
#include <Mt/Separator.h>
#include <iostream.h>
#include <Geometry/Vector.h>
#include <GL/glu.h>
#include <CallbackCloners.h>
#include <Math/MiscMath.h>
#include <Geometry/BBox.h>
#include <stdio.h>
#include <string.h>
#include <X11/keysym.h>

extern MtXEventLoop* evl;

typedef void (Roe::*MouseHandler)(int action, int x, int y,
				  int x_root, int y_root);
#define BUTTON_DOWN 0
#define BUTTON_UP 1
#define BUTTON_MOTION 2
#define SHIFT_MASK 1
#define CONTROL_MASK 2
#define META_MASK 4

struct MouseHandlerData {
    MouseHandler handler;
    clString title;
    MouseHandlerData(MouseHandler, const clString&);
    ~MouseHandlerData();
};

MouseHandlerData::MouseHandlerData(MouseHandler handler, const clString& title)
: handler(handler), title(title)
{
}

MouseHandlerData::~MouseHandlerData()
{
}

static MouseHandlerData mode_translate(&Roe::mouse_translate, "translate");
static MouseHandlerData mode_scale(&Roe::mouse_scale, "scale");
static MouseHandlerData mode_rotate(&Roe::mouse_rotate, "rotate");
static MouseHandlerData mode_pick(&Roe::mouse_pick, "pick");

static MouseHandlerData* mouse_handlers[8][3] = {
    &mode_translate, 	// No modifiers, button 1
    &mode_scale,       	// No modifiers, button 2
    &mode_rotate,	// No modifiers, button 3
    &mode_pick,		// Shift, button 1
    0,			// Shift, button 2
    0,			// Shift, button 3
    0,			// Control, button 1
    0,			// Control, button 2
    0,			// Control, button 3
    0,			// Control+Shift, button 1
    0,			// Control+Shift, button 2
    0,			// Control+Shift, button 3
    0,			// Alt, button 1
    0,			// Alt, button 2
    0,			// Alt, button 3
    0,			// Alt+Shift, button 1
    0,			// Alt+Shift, button 2
    0,			// Alt+Shift, button 3
    0,			// Alt+Control, button 1
    0,			// Alt+Control, button 2
    0,			// Alt+Control, button 3
    0,			// Alt+Control+Shift, button 1
    0,			// Alt+Control+Shift, button 2
    0,			// Alt+Control+Shift, button 3
};

GeomItem::GeomItem() {
}

GeomItem::~GeomItem() {
    delete btn;
}

Roe::Roe(Salmon* s, double *m) {
    haveInheritMat=1;
    for (int i=0; i<16; i++)
	inheritMat[i]=m[i];
    RoeInit(s);
}

Roe::Roe(Salmon* s) {
    haveInheritMat=0;
    RoeInit(s);
}

void Roe::RoeInit(Salmon* s) {
    evl->lock();
    modifier_mask=0;
    doneInit=0;
    old_fh=-1;
    modefont=0;
    buttons_exposed=0;
    manager=s;
    drawinfo=new DrawInfo;
    drawinfo->drawtype=DrawInfo::Gouraud;

    firstGen=False;
    dialog=new DialogShellC;
    dialog->SetAllowShellResize(true);
    dialog->SetDeleteResponse(XmDESTROY);
    new MotifCallback<Roe>FIXCB(dialog, XmNdestroyCallback,
				   &manager->mailbox, this,
				   &Roe::destroyWidgetCB, 0, 0);    
    dialog->SetWidth(600);
    dialog->SetHeight(400);
    dialog->Create("sci", "sci", evl->get_display());

    wholeWin=new RowColumnC;
    wholeWin->SetOrientation(XmHORIZONTAL);
    wholeWin->Create(*dialog, "wholeWin");

    left=new RowColumnC;
    left->SetOrientation(XmVERTICAL);
    left->Create(*wholeWin, "left");

    right=new RowColumnC;
    right->SetOrientation(XmVERTICAL);
    right->Create(*wholeWin, "right");

    gr_frame=new FrameC;
    gr_frame->SetShadowType(XmSHADOW_IN);
    gr_frame->Create(*left, "frame");

    graphics=new GLwMDrawC;
    graphics->SetWidth(400);
    graphics->SetHeight(300);
    graphics->SetRgba(True);
    graphics->SetDoublebuffer(True);
    graphics->SetNavigationType(XmSTICKY_TAB_GROUP);
    graphics->SetTraversalOn(True);
    new MotifCallback<Roe>FIXCB(graphics, GLwNexposeCallback,
				&manager->mailbox, this,
				&Roe::redrawCB,
				0, 0);
    new MotifCallback<Roe>FIXCB(graphics, GLwNginitCallback,
				&manager->mailbox, this,
				&Roe::initCB,
				0, 0);
    new MotifCallback<Roe>FIXCB(graphics, GLwNinputCallback,
				&manager->mailbox, this,
				&Roe::eventCB,
				0, &CallbackCloners::gl_clone);
    graphics->Create(*gr_frame, "opengl_viewer");
    char* translations="<Enter>: glwInput()";
    XtAugmentTranslations(*graphics, XtParseTranslationTable(translations));

    controls=new RowColumnC;
    controls->SetOrientation(XmHORIZONTAL);
    controls->Create(*left, "controls");
    
    objBox=new RowColumnC;
    objBox->SetOrientation(XmVERTICAL);
    objBox->Create(*right, "objBox");
    
    objLabel=new LabelC;
    objLabel->Create(*objBox, "Objects");
    objSep=new SeparatorC;
    objSep->Create(*objBox, "objSep");

    objScroll=new ScrolledWindowC;
    objScroll->SetScrollingPolicy(XmAUTOMATIC);
    objScroll->Create(*right, "objects");

    objRC=new RowColumnC;
    objRC->SetOrientation(XmVERTICAL);
    objRC->Create(*objScroll, "objRC");

    shadeBox=new RowColumnC;
    shadeBox->SetOrientation(XmVERTICAL);
    shadeBox->Create(*right, "shadeBox");
    
    shadeLabel=new LabelC;
    shadeLabel->Create(*shadeBox, "Shading");
    shadeSep=new SeparatorC;
    shadeSep->Create(*shadeBox, "objSep");

    shadeRC=new RowColumnC;
    shadeRC->SetOrientation(XmVERTICAL);
    shadeRC->SetRadioAlwaysOne(True);
    shadeRC->SetRadioBehavior(True);
    shadeRC->Create(*shadeBox, "shadeRC");
    wire=new ToggleButtonC;
    new MotifCallback<Roe>FIXCB(wire, XmNvalueChangedCallback,
				&manager->mailbox, this,
				&Roe::wireCB,
				0, 0);
    wire->Create(*shadeRC, "Wire");
    flat=new ToggleButtonC;
    new MotifCallback<Roe>FIXCB(flat, XmNvalueChangedCallback,
				&manager->mailbox, this,
				&Roe::flatCB,
				0, 0);
    flat->Create(*shadeRC, "Flat");
    gouraud=new ToggleButtonC;
    new MotifCallback<Roe>FIXCB(gouraud, XmNvalueChangedCallback,
				&manager->mailbox, this,
				&Roe::gouraudCB,
				0, 0);
    gouraud->SetSet(True);
    gouraud->Create(*shadeRC, "Gouraud");
    phong=new ToggleButtonC;
    new MotifCallback<Roe>FIXCB(phong, XmNvalueChangedCallback,
			&manager->mailbox, this,
				&Roe::phongCB,
				0, 0);
    phong->Create(*shadeRC, "Phong");

    lightBox=new RowColumnC;
    lightBox->SetOrientation(XmVERTICAL);
    lightBox->Create(*right, "objBox");
    
    lightLabel=new LabelC;
    lightLabel->Create(*lightBox, "Lighting");
    lightSep=new SeparatorC;
    lightSep->Create(*lightBox, "lightSep");

    lightScroll=new ScrolledWindowC;
    lightScroll->SetScrollingPolicy(XmAUTOMATIC);
    lightScroll->Create(*lightBox, "lightScroll");

    lightRC=new RowColumnC;
    lightRC->SetOrientation(XmVERTICAL);
    lightRC->Create(*lightScroll, "lightRC");

    ambient=new ToggleButtonC;
    new MotifCallback<Roe>FIXCB(ambient, XmNvalueChangedCallback,
				&manager->mailbox, this,
				&Roe::ambientCB,
				0, 0);
    ambient->SetSet(True);
    ambient->Create(*lightRC, "Ambient");
    point1=new ToggleButtonC;
    new MotifCallback<Roe>FIXCB(point1, XmNvalueChangedCallback,
				&manager->mailbox, this,
				&Roe::point1CB,
				0, 0);
    point1->SetSet(True);
    point1->Create(*lightRC, "Point1");
    options=new RowColumnC;
    options->SetOrientation(XmHORIZONTAL);
    options->Create(*left, "options");

    viewRC=new RowColumnC;
    viewRC->SetOrientation(XmVERTICAL);
    viewRC->Create(*options, "view");

    autoView=new PushButtonC;
    new MotifCallback<Roe>FIXCB(autoView, XmNactivateCallback,
				&manager->mailbox, this,
				&Roe::autoViewCB,
				0, 0);
    autoView->Create(*viewRC, "Autoview");
    setHome=new PushButtonC;
    new MotifCallback<Roe>FIXCB(setHome, XmNactivateCallback,
				&manager->mailbox, this,
				&Roe::setHomeCB,
				0, 0);
    setHome->Create(*viewRC, "Set Home View");
    goHome=new PushButtonC;
    new MotifCallback<Roe>FIXCB(goHome, XmNactivateCallback,
				&manager->mailbox, this,
				&Roe::goHomeCB,
				0, 0);
    goHome->Create(*viewRC, "Go Home");


    buttons=new DrawingAreaC;
    buttons->SetWidth(200);
    buttons->SetResizePolicy(XmRESIZE_GROW);
    new MotifCallback<Roe>FIXCB(buttons, XmNexposeCallback,
				&manager->mailbox, this,
				&Roe::redraw_buttons,
				0, 0);
    buttons->Create(*options, "buttons");
    gc=XCreateGC(XtDisplay(*buttons), XtWindow(*buttons), 0, 0);
    ColorManager* cm=manager->netedit->color_manager;
    mod_colors[0]=new XQColor(cm, "black");
    mod_colors[1]=new XQColor(cm, "red");
    mod_colors[2]=new XQColor(cm, "purple");
    mod_colors[3]=new XQColor(cm, "blue");
    mod_colors[4]=new XQColor(cm, "green");
    mod_colors[5]=new XQColor(cm, "orange");
    mod_colors[6]=new XQColor(cm, "yellow");
    mod_colors[7]=new XQColor(cm, "white");
    
    spawnRC=new RowColumnC;
    spawnRC->SetOrientation(XmVERTICAL);
    spawnRC->Create(*options, "spawn");

    spawnCh=new PushButtonC;
    new MotifCallback<Roe>FIXCB(spawnCh, XmNactivateCallback,
				&manager->mailbox, this,
				&Roe::spawnChCB,
				0, 0);
    spawnCh->Create(*spawnRC, "Spawn Child");
    spawnInd=new PushButtonC;
    new MotifCallback<Salmon>FIXCB(spawnInd, XmNactivateCallback,
				&manager->mailbox, manager,
				&Salmon::spawnIndCB,
				0, 0);
    spawnInd->Create(*spawnRC, "Spawn Independent");
    
    evl->unlock();
}

void Roe::redrawCB(CallbackData*, void*){
    if(!doneInit)
	initCB(0, 0);
    redrawAll();
}

void Roe::initCB(CallbackData*, void*) {
    XVisualInfo* vi;
    graphics->GetVisualInfo(&vi);
    graphics->GetValues();
    // Create a GLX context
    evl->lock();
    cx = glXCreateContext(XtDisplay(*graphics), vi, 0, GL_TRUE);

    make_current();

    // set the view
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(90, 1.33, 1, 10);
    glMatrixMode(GL_MODELVIEW);
    if (haveInheritMat) {
	glLoadMatrixd(inheritMat);
    } else {
	glLoadIdentity();
	gluLookAt(2,2,5,2,2,2,0,1,0);
	glGetDoublev(GL_MODELVIEW_MATRIX, inheritMat);
    }

    GLfloat light_position[] = { 3,3,-100,1};
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);
    evl->unlock();
    doneInit=1;
}

void Roe::make_current() {
    evl->lock();
    GLwDrawingAreaMakeCurrent(*graphics, cx);
    evl->unlock();
}

void Roe::itemAdded(GeomObj *g, char *name) {
    GeomItem *item;
    item= new GeomItem;
    ToggleButtonC *bttn;
    bttn = new ToggleButtonC;
    item->btn=bttn;
    item->vis=1;
    item->geom=g;
    strcpy(item->name, name);
    geomItemA.add(item);
    new MotifCallback<Roe>FIXCB(item->btn, XmNvalueChangedCallback,
				&manager->mailbox, this,
				&Roe::itemCB,
				(void *) item, 0);
    item->btn->SetSet(True);
    item->btn->Create(*objRC, name);
    for (int i=0; i<kids.size(); i++) {
	kids[i]->itemAdded(g, name);
    }
}

void Roe::itemDeleted(GeomObj *g) {
    for (int i=0; i<geomItemA.size(); i++) {
	if (geomItemA[i]->geom == g) {
	    delete (geomItemA[i]->btn);
	    geomItemA.remove(i);
	}
    }
    for (i=0; i<kids.size(); i++) {
	kids[i]->itemDeleted(g);
    }
}


void Roe::redrawAll()
{
    if (doneInit) {
	// clear screen
	evl->lock();
        make_current();  
	glClearColor(0,0,0,1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	drawinfo->push_matl(manager->default_matl);
	HashTableIter<int,HashTable<int, GeomObj*>*> iter(&manager->portHash);
	for (iter.first(); iter.ok(); ++iter) {
	    HashTable<int, GeomObj*>* serHash=iter.get_data();
	    HashTableIter<int, GeomObj*> serIter(serHash);
	    for (serIter.first(); serIter.ok(); ++serIter) {
		GeomObj *geom=serIter.get_data();
		for (int i=0; i<geomItemA.size(); i++)
		    if (geomItemA[i]->geom == geom)
			if (geomItemA[i]->vis)
			    geom->draw(drawinfo);
	    }	
	}
	drawinfo->pop_matl();
	GLwDrawingAreaSwapBuffers(*graphics);
	for (int i=0; i<kids.size(); i++) {
	    kids[i]->redrawAll();
	}
	evl->unlock();       
    }
}

void Roe::printLevel(int level, int&flag) {
    if (level == 0) {
	flag=1;
	cerr << "* ";
    } else {
	for (int i=0; i<kids.size(); i++) {
	    kids[i]->printLevel(level-1, flag);
	}
    }
}
 
// need to fill this in!   
void Roe::itemCB(CallbackData*, void *gI) {
    GeomItem *g = (GeomItem *)gI;
    for (int i=0; i<geomItemA.size(); i++) {
	if (geomItemA[i]->geom == g->geom) {
	    if (geomItemA[i]->vis) {
		geomItemA[i]->vis=0;
	    } else {
		geomItemA[i]->vis=1;
	    }
	}
    }
    redrawAll();
}

void Roe::destroyWidgetCB(CallbackData*, void*)
{
    // can't close the only window -- this doesn't seem to work, though...
    if (firstGen && (manager->topRoe.size() == 1) && (kids.size()==0)) 
	return;
    else
	delete this;
}

void Roe::spawnChCB(CallbackData*, void*)
{
  double mat[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, mat);

/*  for (int i=0;i<16;i++)
      cerr << mat[i] << " ";
  cerr << "\n";
*/
  kids.add(new Roe(manager, mat));
  kids[kids.size()-1]->SetParent(this);
  for (int i=0; i<geomItemA.size(); i++)
      kids[kids.size()-1]->itemAdded(geomItemA[i]->geom, geomItemA[i]->name);
  kids[kids.size()-1]->redrawAll();
//  manager->printFamilyTree();

}
    
Roe::~Roe()
{
    delete dialog;
    delete wholeWin;
    delete left;
    delete right;
    delete graphics;
    delete controls;
    delete objBox;
    delete objLabel;
    delete objSep;
    delete objScroll;
    delete objRC;
    delete shadeBox;
    delete shadeLabel;
    delete shadeSep;
    delete shadeRC;
    delete wire;
    delete flat;
    delete phong;
    delete gouraud;
    delete lightBox;
    delete lightLabel;
    delete lightScroll;
    delete lightSep;
    delete lightRC;
    delete ambient;
    delete point1;
    delete options;
    delete viewRC;
    delete autoView;
    delete setHome;
    delete goHome;
    delete buttons;
    delete spawnRC;
    delete spawnCh;
    delete spawnInd;
    delete form;
    delete gr_frame;

    delete drawinfo;

    for (int i=0; i<geomItemA.size(); i++)
	delete geomItemA[i];
    geomItemA.remove_all();

    // tell my parent to delete me from their kid list
    if (firstGen) {
	manager->delTopRoe(this);
    } else {
	parent->deleteChild(this);
    }

    // now take care of the kids!  If I'm first generation, add them
    // to the Salmon's topRoe and delete myself from the Salmon's topRoe
    // Otherwise, give them to my parents and delete myself from my parents
    // Don't forget to set their firstGen, and parent variables accordingly
    if (firstGen) {
	for (int i=0; i<kids.size(); i++) {
	    kids[i]->SetTop();
	    manager->addTopRoe(kids[i]);
	}
    } else {
	for (int i=0; i<kids.size(); i++) {
	    kids[i]->SetParent(parent);
	    parent->addChild(kids[i]);
	}
    }
//    manager->printFamilyTree();
}

void Roe::SetParent(Roe *r)
{
  parent = r;
}

void Roe::SetTop()
{
  firstGen=True;
}

void Roe::addChild(Roe *r)
{
    kids.add(r);
}

// self-called method
void Roe::deleteChild(Roe *r)
{
    for (int i=0; i<kids.size(); i++)
	if (r==kids[i]) kids.remove(i);
}
void Roe::wireCB(CallbackData*, void*)
{
    drawinfo->drawtype=DrawInfo::WireFrame;
    drawinfo->current_matl=0;
    make_current();
    glDisable(GL_LIGHTING);
    redrawAll();
}

void Roe::flatCB(CallbackData*, void*)
{
    drawinfo->drawtype=DrawInfo::Flat;
    drawinfo->current_matl=0;
    make_current();
    glDisable(GL_LIGHTING);
    redrawAll();
}

void Roe::gouraudCB(CallbackData*, void*)
{
    drawinfo->drawtype=DrawInfo::Gouraud;
    drawinfo->current_matl=0;
    make_current();
    glEnable(GL_LIGHTING);
    redrawAll();
}

void Roe::phongCB(CallbackData*, void*)
{
    drawinfo->drawtype=DrawInfo::Phong;
    drawinfo->current_matl=0;
    make_current();
    glEnable(GL_LIGHTING);
    redrawAll();
}

void Roe::ambientCB(CallbackData*, void*)
{
    NOT_FINISHED("Roe::ambientCB");
}
void Roe::point1CB(CallbackData*, void*)
{
    make_current();
    if (!glIsEnabled(GL_LIGHT0)) {
	glEnable(GL_LIGHT0);
    } else {
	glDisable(GL_LIGHT0);
    }
    redrawAll();
}
void Roe::goHomeCB(CallbackData*, void*)
{
    make_current();
    glLoadMatrixd(inheritMat);
    redrawAll();
}
void Roe::autoViewCB(CallbackData*, void*)
{
    BBox bbox;
    int ok=0;
    HashTableIter<int,HashTable<int, GeomObj*>*> iter(&manager->portHash);
    for (iter.first(); iter.ok(); ++iter) {
	HashTable<int, GeomObj*>* serHash=iter.get_data();
	HashTableIter<int, GeomObj*> serIter(serHash);
	for (serIter.first(); serIter.ok(); ++serIter) {
	    GeomObj *geom=serIter.get_data();
	    for (int i=0; i<geomItemA.size(); i++)
		if (geomItemA[i]->geom == geom)
		    if (geomItemA[i]->vis) {
			bbox.extend(geom->bbox());
			ok=1;
		    }
	}		
    }	
    if (!ok) return;
    Point lookat(bbox.center());
    lookat.z(bbox.max().z());
    double xwidth=lookat.x()-bbox.min().x();
    double ywidth=lookat.y()-bbox.min().y();
    double dist=Max(xwidth, ywidth);
    make_current();
    evl->lock();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(90, 1.33, 1, 1000);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(lookat.x(), lookat.y(), lookat.z()+dist, lookat.x(), lookat.y(), lookat.z(), 0, 1, 0);
    evl->unlock();
    redrawAll();
}    

void Roe::setHomeCB(CallbackData*, void*)
{
    make_current();
    glGetDoublev(GL_MODELVIEW_MATRIX, inheritMat);
}

Roe::Roe(const Roe& copy)
{
    NOT_FINISHED("Roe::Roe");
}

void Roe::rotate(double angle, Vector v)
{
    make_current();
    glRotated(angle,v.x(), v.y(), v.z());
}

void Roe::translate(Vector v)
{
    make_current();
    glTranslated(v.x(), v.y(), v.z());
}

void Roe::scale(Vector v)
{
    make_current();
    glScaled(v.x(), v.y(), v.z());
}

void Roe::eventCB(CallbackData* cbdata, void*)
{
    XEvent* event=cbdata->get_event();
    switch(event->type){
    case EnterNotify:
	evl->lock();
	XmProcessTraversal(*graphics, XmTRAVERSE_CURRENT);
	evl->unlock();
	return;
    case KeyPress:
	if(event->xkey.state & (Button1Mask|Button2Mask|Button3Mask)){
	    // Skip it...
	} else{
	    int mask=0;
	    if(event->xkey.state & ControlMask)
		mask|=CONTROL_MASK;
	    if(event->xkey.state & ShiftMask)
		mask|=SHIFT_MASK;
	    if(event->xkey.state & (Mod1Mask|Mod2Mask|Mod3Mask|Mod4Mask|Mod5Mask))
		mask|=META_MASK;
	    switch(XLookupKeysym(&event->xkey, 0)){
	    case XK_Shift_L:
	    case XK_Shift_R:
		mask|=SHIFT_MASK;
		break;
	    case XK_Control_L:
	    case XK_Control_R:
		mask|=CONTROL_MASK;
		break;
	    case XK_Meta_L:
	    case XK_Meta_R:
	    case XK_Alt_L:
	    case XK_Alt_R:
		mask|=META_MASK;
		break;
	    }
	    
	    if(mask != modifier_mask){
		modifier_mask=mask;
		update_modifier_widget();
	    }
	}
	break;
    case KeyRelease:
	if(event->xkey.state & (Button1Mask|Button2Mask|Button3Mask)){
	    // Skip it...
	} else{
	    int mask=0;
	    if(event->xkey.state & ControlMask)
		mask|=CONTROL_MASK;
	    if(event->xkey.state & ShiftMask)
		mask|=SHIFT_MASK;
	    if(event->xkey.state & (Mod1Mask|Mod2Mask|Mod3Mask|Mod4Mask|Mod5Mask))
		mask|=META_MASK;
	    switch(XLookupKeysym(&event->xkey, 0)){
	    case XK_Shift_L:
	    case XK_Shift_R:
		mask&=~SHIFT_MASK;
		break;
	    case XK_Control_L:
	    case XK_Control_R:
		mask&=~CONTROL_MASK;
		break;
	    case XK_Meta_L:
	    case XK_Meta_R:
	    case XK_Alt_L:
	    case XK_Alt_R:
		mask&=~META_MASK;
		break;
	    }
	    if(mask != modifier_mask){
		modifier_mask=mask;
		update_modifier_widget();
	    }
	}
	break;
    case ButtonPress:
	{
	    switch(event->xbutton.button){
	    case Button1:
		last_btn=1;
		break;
	    case Button2:
		last_btn=2;
		break;
	    case Button3:
		last_btn=3;
		break;
	    default:
		last_btn=1;
		break;
	    }
	    int mask=0;
	    if(event->xbutton.state & ControlMask)
		mask|=CONTROL_MASK;
	    if(event->xbutton.state & ShiftMask)
		mask|=SHIFT_MASK;
	    if(event->xbutton.state & (Mod1Mask|Mod2Mask|Mod3Mask|Mod4Mask|Mod5Mask))
		mask|=META_MASK;
	    if(mask != modifier_mask){
		modifier_mask=mask;
		update_modifier_widget();
	    }
	    MouseHandler handler=mouse_handlers[modifier_mask][last_btn-1]->handler;
	    if(handler){
		(this->*handler)(BUTTON_DOWN,
				 event->xbutton.x, event->xbutton.y,
				 event->xbutton.x_root, event->xbutton.y_root);
	    }
	}
	break;
    case ButtonRelease:
	{
	    switch(event->xbutton.button){
	    case Button1:
		last_btn=1;
		break;
	    case Button2:
		last_btn=2;
		break;
	    case Button3:
		last_btn=3;
		break;
	    default:
		last_btn=1;
		break;
	    }
	    MouseHandler handler=mouse_handlers[modifier_mask][last_btn-1]->handler;
	    if(handler){
		(this->*handler)(BUTTON_UP,
				 event->xbutton.x, event->xbutton.y,
				 event->xbutton.x_root, event->xbutton.y_root);
	    }
	    int mask=0;
	    if(event->xbutton.state & ControlMask)
		mask|=CONTROL_MASK;
	    if(event->xbutton.state & ShiftMask)
		mask|=SHIFT_MASK;
	    if(event->xbutton.state & (Mod1Mask|Mod2Mask|Mod3Mask|Mod4Mask|Mod5Mask))
		mask|=META_MASK;
	    if(mask != modifier_mask){
		modifier_mask=mask;
		update_modifier_widget();
	    }
	}
	break;
    case MotionNotify:
	{
	    MouseHandler handler=mouse_handlers[modifier_mask][last_btn-1]->handler;
	    if(handler){
		(this->*handler)(BUTTON_MOTION,
				 event->xmotion.x, event->xmotion.y,
				 event->xmotion.x_root, event->xmotion.y_root);
	    }
	}
	break;
    default:
	cerr << "Unknown event..\n";
	break;
    }
}

void Roe::mouse_translate(int action, int x, int y, int, int)
{
    switch(action){
    case BUTTON_DOWN:
	last_x=x;
	last_y=y;
	update_mode_string("translate: ");
	break;
    case BUTTON_MOTION:
	{
	    double xmtn=last_x-x;
	    double ymtn=last_y-y;
	    xmtn/=10;
	    ymtn/=10;
	    last_x = x;
	    last_y = y;
	    make_current();
	    glTranslated(-xmtn, ymtn, 0);
	    for (int i=0; i<kids.size(); i++)
		kids[i]->translate(Vector(-xmtn, ymtn, 0));
	    redrawAll();
	    update_mode_string(clString("translate: ")+to_string(xmtn)
			       +", "+to_string(ymtn));
	}
	break;
    case BUTTON_UP:
	update_mode_string("");
	break;
    }
}

void Roe::mouse_scale(int action, int x, int y, int, int)
{
    switch(action){
    case BUTTON_DOWN:
	last_x=x;
	last_y=y;
	break;
    case BUTTON_MOTION:
	{
	    double scl;
	    double xmtn=last_x-x;
	    double ymtn=last_y-y;
	    xmtn/=30;
	    ymtn/=30;
	    last_x = x;
	    last_y = y;
	    make_current();
	    if (Abs(xmtn)>Abs(ymtn)) scl=xmtn; else scl=ymtn;
	    glScaled(1+scl, 1+scl, 1+scl);
	    for (int i=0; i<kids.size(); i++)
		kids[i]->scale(Vector(1+scl, 1+scl, 1+scl));
	    redrawAll();
	}
	break;
    }
}

void Roe::mouse_rotate(int action, int x, int y, int, int)
{
    switch(action){
    case BUTTON_DOWN:
	{
	    last_x=x;
	    last_y=y;
	    bb.reset();
	    HashTableIter<int,HashTable<int, GeomObj*>*> iter(&manager->portHash);
	    for (iter.first(); iter.ok(); ++iter) {
		HashTable<int, GeomObj*>* serHash=iter.get_data();
		HashTableIter<int, GeomObj*> serIter(serHash);
		for (serIter.first(); serIter.ok(); ++serIter) {
		    GeomObj *geom=serIter.get_data();
		    for (int i=0; i<geomItemA.size(); i++)
			if (geomItemA[i]->geom == geom)
			    if (geomItemA[i]->vis)
				bb.extend(geom->bbox());
		}		
	    }	
	}
	break;
    case BUTTON_MOTION:
	{
	    double xmtn=last_x-x;
	    double ymtn=last_y-y;
	    last_x = x;
	    last_y = y;
	    make_current();
	    Point cntr(bb.center());
	    glTranslated(cntr.x(), cntr.y(), cntr.z());
	    glRotated(xmtn*xmtn+ymtn*ymtn,-ymtn,-xmtn,0);
	    glTranslated(-cntr.x(), -cntr.y(), -cntr.z());
	    for (int i=0; i<kids.size(); i++) {
		kids[i]->translate(Vector(cntr.x(), cntr.y(), cntr.z()));
		kids[i]->rotate(xmtn*xmtn+ymtn*ymtn,Vector(-ymtn,-xmtn,0));
		kids[i]->translate(Vector(-cntr.x(), -cntr.y(), -cntr.z()));
	    }
	    redrawAll();
	}
	break;
    case BUTTON_UP:
	break;
    }
}

void Roe::mouse_pick(int action, int x, int y, int, int)
{
    NOT_FINISHED("mouse_pick");
}

void Roe::update_mode_string(const clString& ms)
{
    mode_string=ms;
    update_modifier_widget();
}

void Roe::update_modifier_widget()
{
    evl->lock();
    if(!buttons_exposed)return;
    Window w=XtWindow(*buttons);
    Display* dpy=XtDisplay(*buttons);
    XClearWindow(dpy, w);
    XSetForeground(dpy, gc, mod_colors[modifier_mask]->pixel());
    Dimension h;
    buttons->GetHeight(&h);
    buttons->GetValues();
    int fh=h/5;
    if(fh != old_fh || modefont==0){
	if(modefont)XUnloadFont(dpy, modefont);
	old_fh=fh;
	char pattern[1000];
	sprintf(pattern, "screen14");
	int acount;
	char** fontnames=XListFonts(dpy, pattern, 1, &acount);
	if(acount==1){
	    modefont=XLoadFont(dpy, fontnames[0]);
	    XFreeFontNames(fontnames);
	    XSetFont(dpy, gc, modefont);
	}
    }
    int fh2=fh/2;
    int fh4=fh2/2;
    int wid=fh2+fh4;
    XFillArc(dpy, w, gc, fh2,       fh2, fh2, fh2, 0, 180*64);
    XFillArc(dpy, w, gc, fh2+wid,   fh2, fh2, fh2, 0, 180*64);
    XFillArc(dpy, w, gc, fh2+2*wid, fh2, fh2, fh2, 0, 180*64);
    XFillArc(dpy, w, gc, fh2,       fh,  fh2, fh2, 180*64, 180*64);
    XFillArc(dpy, w, gc, fh2+wid,   fh,  fh2, fh2, 180*64, 180*64);
    XFillArc(dpy, w, gc, fh2+2*wid, fh,  fh2, fh2, 180*64, 180*64);
    XFillRectangle(dpy, w, gc, fh2,       fh2+fh4, fh2+1, fh2+2);
    XFillRectangle(dpy, w, gc, fh2+wid,   fh2+fh4, fh2+1, fh2+2);
    XFillRectangle(dpy, w, gc, fh2+2*wid, fh2+fh4, fh2+1, fh2+2);

    int toff=wid*3+fh;
    XDrawLine(dpy, w, gc, fh2+fh4, fh+fh4, fh2+fh4, 2*fh);
    XDrawLine(dpy, w, gc, fh2+fh4, 2*fh  , toff-fh4, 2*fh);
    XDrawLine(dpy, w, gc, fh2+wid+fh4, fh+fh4, fh2+wid+fh4, 3*fh);
    XDrawLine(dpy, w, gc, fh2+wid+fh4, 3*fh,   toff-fh4, 3*fh);
    XDrawLine(dpy, w, gc, fh2+2*wid+fh4, fh+fh4, fh2+2*wid+fh4, 4*fh);
    XDrawLine(dpy, w, gc, fh2+2*wid+fh4, 4*fh, toff-fh4, 4*fh);

    XSetForeground(dpy, gc, BlackPixelOfScreen(XtScreen(*buttons)));
    XDrawString(dpy, w, gc, toff, fh+fh2, mode_string(), mode_string.len());
    clString b1_string(mouse_handlers[modifier_mask][0]?
		       mouse_handlers[modifier_mask][0]->title
		       :clString(""));
    clString b2_string(mouse_handlers[modifier_mask][1]?
		       mouse_handlers[modifier_mask][1]->title
		       :clString(""));
    clString b3_string(mouse_handlers[modifier_mask][2]?
		       mouse_handlers[modifier_mask][2]->title
		       :clString(""));
    XDrawString(dpy, w, gc, toff, 2*fh+fh2, b1_string(), b1_string.len());
    XDrawString(dpy, w, gc, toff, 3*fh+fh2, b2_string(), b2_string.len());
    XDrawString(dpy, w, gc, toff, 4*fh+fh2, b3_string(), b3_string.len());
    evl->unlock();
}

void Roe::redraw_buttons(CallbackData*, void*)
{
    buttons_exposed=1;
    update_modifier_widget();
}

