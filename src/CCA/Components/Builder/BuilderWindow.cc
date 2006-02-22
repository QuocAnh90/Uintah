/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2004 Scientific Computing and Imaging Institute,
   University of Utah.

   License for the specific language governing rights and limitations under
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/

#include <wx/event.h>
#include <wx/frame.h>
#include <wx/image.h>
#include <wx/textctrl.h>
#include <wx/laywin.h>
#include <wx/string.h>

#include <vector>
#include <map>

#include <Core/Thread/Thread.h>
#include <Core/Containers/StringUtil.h>

#include <SCIRun/TypeMap.h>
#include <SCIRun/CCA/CCAException.h>

#include <CCA/Components/Builder/BuilderWindow.h>
#include <CCA/Components/Builder/MiniCanvas.h>
#include <CCA/Components/Builder/NetworkCanvas.h>
#include <CCA/Components/Builder/ComponentIcon.h>

namespace GUIBuilder {

using namespace SCIRun;

MenuTree::MenuTree(BuilderWindow* builder, const std::string &url) : builder(builder), url(url), id(0)
{
}

MenuTree::~MenuTree()
{
  for (std::map<std::string, MenuTree*>::iterator iter = child.begin(); iter != child.end(); iter++) {
    delete iter->second;
  }
  if (! cd.isNull()) {
    // event handler cleanup
    this->Disconnect(id,  wxEVT_COMMAND_MENU_SELECTED,
		     wxCommandEventHandler(MenuTree::OnInstantiateComponent));
    builder->RemoveEventHandler(this);
  }
}

void MenuTree::add(const std::vector<std::string>& name, int nameindex,
           const sci::cca::ComponentClassDescription::pointer& desc,
           const std::string& fullname)
{
  if (nameindex == (int) name.size()) {
    if ( !cd.isNull() ) {
      // warning - should be displayed?
      std::cerr << "Duplicate component: " << fullname << '\n';
    } else {
      cd = desc;
      id = BuilderWindow::GetNextID();
    }
  } else {
    const std::string& n = name[nameindex];
    std::map<std::string, MenuTree*>::iterator iter = child.find(n);
    if(iter == child.end()) {
      child[n] = new MenuTree(builder, url);
    }
    child[n]->add(name, nameindex + 1, desc, fullname);
  }
}

// Consolidate component class names from the bottom up.
void MenuTree::coalesce()
{
  for (std::map<std::string, MenuTree*>::iterator iter = child.begin();
       iter != child.end(); iter++) {
    MenuTree* c = iter->second;
    while (c->child.size() == 1) {
      std::map<std::string, MenuTree*>::iterator grandchild = c->child.begin();
      std::string newname = iter->first + "." + grandchild->first;

      MenuTree* gc = grandchild->second;
      c->child.clear(); // So that grandchild won't get deleted...
      delete c;

      child.erase(iter);
      child[newname] = gc;
      iter = child.begin();
      c = gc;
    }
    c->coalesce();
  }
}

void MenuTree::populateMenu(wxMenu* menu)
{
  for (std::map<std::string, MenuTree*>::iterator iter = child.begin();
       iter != child.end(); iter++) {
    if (iter->second->cd.isNull()) {
      wxMenu* submenu = new wxMenu(wxT(""), wxMENU_TEAROFF);
      //submenu->setFont(builder->font());
      iter->second->populateMenu(submenu);
      menu->Append(ID_MENU_COMPONENTS, wxT(iter->first), submenu);
    } else {
      builder->PushEventHandler(iter->second);
      menu->Append(iter->second->id, wxT(iter->first), wxT(iter->first));
      iter->second->Connect(iter->second->id, wxEVT_COMMAND_MENU_SELECTED,
			    wxCommandEventHandler(MenuTree::OnInstantiateComponent));
    }
  }
}

void MenuTree::clear()
{
  child.clear();
}

void MenuTree::OnInstantiateComponent(wxCommandEvent& event)
{
  // this shouldn't happen
  if (cd.isNull()) {
    // error should be logged
    std::cerr << "Error: null component description!" << std::endl;
    return;
  }
  builder->InstantiateComponent(cd);
}



const wxColor BuilderWindow::BACKGROUND_COLOUR(0, 51, 102);

// Event table
BEGIN_EVENT_TABLE(BuilderWindow, wxFrame)
  EVT_MENU(wxID_ABOUT, BuilderWindow::OnAbout)
  EVT_MENU(wxID_EXIT, BuilderWindow::OnQuit)
  EVT_MENU(ID_MENU_TEST, BuilderWindow::OnTest)
  //EVT_MENU(MenuTree::ID_MENU_COMPONENTS, MenuTree::OnInstantiateComponent)
  EVT_SIZE(BuilderWindow::OnSize)
  EVT_SASH_DRAGGED_RANGE(ID_WINDOW_LEFT, ID_WINDOW_BOTTOM, BuilderWindow::OnSashDrag)
END_EVENT_TABLE()

IMPLEMENT_DYNAMIC_CLASS(BuilderWindow, wxFrame)

int BuilderWindow::IdCounter = BuilderWindow::ID_BUILDERWINDOW_HIGHEST;

BuilderWindow::BuilderWindow(const sci::cca::BuilderComponent::pointer& bc, wxWindow *parent) : builder(bc)
{
std::cerr << "BuilderWindow::BuilderWindow(..): from thread " << Thread::self()->getThreadName() << " in framework " << builder->getFrameworkURL() << std::endl;

  Init();
  Create(parent, wxID_ANY);
}

bool BuilderWindow::Create(wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style, const wxString& name)
{
  if (!wxFrame::Create(parent, id, title, pos, size, style, name)) {
    return false;
  }

  bottomWindow = new wxSashLayoutWindow(this, ID_WINDOW_BOTTOM, wxPoint(0, TOP_HEIGHT), wxSize(WIDTH, BOTTOM_HEIGHT), wxNO_BORDER|wxSW_3D|wxCLIP_CHILDREN);
  bottomWindow->SetDefaultSize(wxSize(MIN, BOTTOM_HEIGHT));
  bottomWindow->SetOrientation(wxLAYOUT_HORIZONTAL);
  bottomWindow->SetAlignment(wxLAYOUT_BOTTOM);
  bottomWindow->SetSashVisible(wxSASH_TOP, true);

  networkCanvas = new NetworkCanvas(builder, this, bottomWindow, ID_NET_WINDOW, wxPoint(0, TOP_HEIGHT), wxSize(WIDTH, BOTTOM_HEIGHT));

  // use wxCLIP_CHILDREN to eliminate flicker on Windows
  // A window to the left of the client window
  leftWindow = new wxSashLayoutWindow(this, ID_WINDOW_LEFT, wxPoint(0, 0), wxSize(MINI_WIDTH, TOP_HEIGHT), wxNO_BORDER|wxSW_3D|wxCLIP_CHILDREN);
  leftWindow->SetDefaultSize(wxSize(MINI_WIDTH, MIN));
  leftWindow->SetOrientation(wxLAYOUT_VERTICAL);
  leftWindow->SetAlignment(wxLAYOUT_LEFT);
  leftWindow->SetBackgroundColour(BACKGROUND_COLOUR);
  leftWindow->SetSashVisible(wxSASH_RIGHT, true);

  // add mini-canvas (scrolled window) to leftWindow
  miniCanvas = new MiniCanvas(leftWindow, networkCanvas, ID_MINI_WINDOW, wxPoint(0, 0), wxSize(MINI_WIDTH, TOP_HEIGHT));

  // A window to the left of the client window
  rightWindow = new wxSashLayoutWindow(this, ID_WINDOW_RIGHT, wxPoint(MINI_WIDTH, 0), wxSize(TEXT_WIDTH, TOP_HEIGHT), wxNO_BORDER|wxSW_3D|wxCLIP_CHILDREN);
  rightWindow->SetDefaultSize(wxSize(TEXT_WIDTH, MIN));
  rightWindow->SetOrientation(wxLAYOUT_VERTICAL);
  rightWindow->SetAlignment(wxLAYOUT_LEFT);
  rightWindow->SetSashVisible(wxSASH_LEFT, false); // resizing in terms of leftWindow only

  url = builder->getFrameworkURL();
  textCtrl = new wxTextCtrl(rightWindow, ID_TEXT_WINDOW, wxT(""), wxPoint(MINI_WIDTH, 0), wxSize(TEXT_WIDTH, TOP_HEIGHT), wxSUNKEN_BORDER|wxTE_MULTILINE|wxTE_READONLY|wxTE_AUTO_URL|wxTE_CHARWRAP);
  (*textCtrl) << "SCIRun2 v 0.1\n";
  (*textCtrl) << "Framework URL: " << url.c_str() << "\n";
  (*textCtrl) << "--------------------\n" << "\n";

  // The "About" item should be in the help menu
  wxMenu *helpMenu = new wxMenu();
  helpMenu->Append(wxID_ABOUT, wxT("&About...\tF1"), wxT("Show about dialog"));

  wxMenu* compWizardMenu = new wxMenu(wxT(""), wxMENU_TEAROFF);
  compWizardMenu->Append(ID_MENU_COMPONENT_WIZARD, wxT("Component Wizard"), wxT("Create component skeleton"));

  wxMenu* fileMenu = new wxMenu();
  fileMenu->Append(ID_MENU_TEST, wxT("&Test\tAlt-T"), wxT("Test component build"));
  fileMenu->AppendSeparator();
  fileMenu->Append(ID_MENU_LOAD, wxT("&Load\tAlt-L"), wxT("Load network file"));
  fileMenu->Append(ID_MENU_INSERT, wxT("&Insert\tAlt-L"), wxT("Insert network file"));
  fileMenu->Append(wxID_SAVE, wxT("&Save\tAlt-S"), wxT("Save network to file"));
  fileMenu->Append(wxID_SAVEAS, wxT("&Save As\tAlt-S"), wxT("Save network to file"));
  fileMenu->AppendSeparator();
  fileMenu->Append(ID_MENU_CLEAR, wxT("&Clear\tAlt-C"), wxT("Clear All"));
  fileMenu->Append(wxID_SELECTALL, wxT("Select &All\tCtrl-A"), wxT("Select All"));
  //fileMenu->Append(ID_MENU_EXECALL, wxT("&Execute All\tCtrl-A"), wxT("Execute All"));
  fileMenu->AppendSeparator();
  fileMenu->Append(ID_MENU_WIZARDS, wxT("&Wizards\tAlt-W"), compWizardMenu);
  fileMenu->AppendSeparator();
  //fileMenu->Append(ID_MENU_ADDINFO, wxT("&Add Info\tAlt-A"), wxT("Add information to?"));
  //fileMenu->AppendSeparator();
  fileMenu->Append(wxID_EXIT, wxT("E&xit\tAlt-X"), wxT("Quit this program"));

  menuBar = new wxMenuBar();
  menuBar->Append(fileMenu, wxT("&File"));

  buildPackageMenus();

  menuBar->Append(helpMenu, wxT("&Help"));
  SetMenuBar(menuBar);

  //SetFont(wxFont(11, wxDEFAULT, wxNORMAL, wxNORMAL, 0, wxT("Sans")));

  statusBar = CreateStatusBar(2, wxST_SIZEGRIP);
  int statusBarWidths[] = { 350, -1 };
  statusBar->SetStatusWidths(2, statusBarWidths);
//     const wxString statusBar_fields[] = {
//         wxT(""),
//         wxT("")
//     };
//     for(int i = 0; i < statusBar->GetFieldsCount(); ++i) {
//         statusBar->SetStatusText(statusBar_fields[i], i);
//     }

  statusBar->SetStatusText("SCIRun2 started", 0);

  return true;
}

void BuilderWindow::Init()
{
std::cerr << "BuilderWindow::Init()" << std::endl;
  //SetFont(wxFont(11, wxDEFAULT, wxNORMAL, wxNORMAL, 0, wxT("Sans")));
  //SetToolTip(wxT("\"Test tooltip\""));
}

bool BuilderWindow::SetBuilder(const sci::cca::BuilderComponent::pointer& bc)
{
  if (builder.isNull()) {
    builder = bc;
    return true;
  }

  return false;
}

BuilderWindow::~BuilderWindow()
{
  std::cerr << "BuilderWindow::~BuilderWindow()" << std::endl;
  // framework shutdown instead!!!
  Thread::exitAll(0);
}

void BuilderWindow::OnAbout(wxCommandEvent &event)
{
  wxString msg;
  msg.Printf(wxT("Hello and welcome to %s"),
             wxVERSION_STRING);

  // show license
  wxMessageBox(msg, wxT("About SCIRun2"),
               wxOK | wxICON_INFORMATION, this);
}

void BuilderWindow::OnQuit(wxCommandEvent &event)
{
  std::cerr << "BuilderWindow::OnQuit(..)" << std::endl;
  // Destroy the frame
  Close();
}

void BuilderWindow::OnSashDrag(wxSashEvent& event)
{
  if (event.GetDragStatus() == wxSASH_STATUS_OUT_OF_RANGE) {
    return;
  }

  switch (event.GetId()) {
  // If both the left and right windows attempt to handle dragging events,
  // the delivery of events can be error-prone, resulting in resizing
  // errors.
  case ID_WINDOW_LEFT:
    leftWindow->SetDefaultSize(wxSize(event.GetDragRect().width, MIN));
    break;
  case ID_WINDOW_BOTTOM:
    bottomWindow->SetDefaultSize(wxSize(MIN, event.GetDragRect().height));
    break;
  }

  wxLayoutAlgorithm layout;
  layout.LayoutFrame(this);

  // Leaves bits of itself behind sometimes
  Refresh();
}

void BuilderWindow::OnSize(wxSizeEvent& WXUNUSED(event))
{
  // recalc. window sizes so that bottom window gets more
  wxLayoutAlgorithm layout;
  layout.LayoutFrame(this);
  Refresh();
}

void BuilderWindow::RedrawMiniCanvas()
{
  miniCanvas->Refresh();
}

void BuilderWindow::OnTest(wxCommandEvent&/* event */)
{
  statusBar->SetStatusText("Build components", 0);
  try {
    sci::cca::ComponentID::pointer helloCid = builder->createInstance("SCIRun.Hello", sci::cca::TypeMap::pointer(0));
    if (! helloCid.isNull()) {
      std::cerr << "wx: Got hello: " << helloCid->getInstanceName() << std::endl;
      networkCanvas->AddIcon(helloCid);
    }

//     sci::cca::ComponentID::pointer helloCid1 = builder->createInstance("SCIRun.Hello");
//     std::cerr << "wx: Got hello: " << helloCid1->getInstanceName() << std::endl;

//     sci::cca::ComponentID::pointer helloCid2 = builder->createInstance("SCIRun.Hello");
//     std::cerr << "wx: Got hello: " << helloCid2->getInstanceName() << std::endl;

    sci::cca::ComponentID::pointer worldCid = builder->createInstance("SCIRun.World", sci::cca::TypeMap::pointer(0));
    if (! worldCid.isNull()) {
      std::cerr << "wx: Got world: " << worldCid->getInstanceName() << std::endl;
      networkCanvas->AddIcon(worldCid);
    }

    sci::cca::ComponentID::pointer pdeDriverCid = builder->createInstance("SCIRun.PDEdriver", sci::cca::TypeMap::pointer(0));
    if (! pdeDriverCid.isNull()) {
      std::cerr << "wx: Got pdeDriver: " << pdeDriverCid->getInstanceName() << std::endl;
      networkCanvas->AddIcon(pdeDriverCid);
    }

//     sci::cca::ComponentID::pointer builderCid = builder->createInstance("SCIRun.Builder");
//     std::cerr << "wx: Got builder: " << builderCid->getInstanceName() << std::endl;

    //sci::cca::ConnectionID::pointer connID = bs->connect(worldCid, "stringport", helloCid, "stringport");
    //Refresh();

    //services->registerUsesPort("goPort", "sci.cca.ports.GoPort", sci::cca::TypeMap::pointer(0));
    //bs->connect(services->getComponentID(), "goPort", helloCid, "go");

    //sci::cca::Port::pointer p = services->getPort("goPort");
    //sci::cca::ports::GoPort::pointer goPort = pidl_cast<sci::cca::ports::GoPort::pointer>(p);
    //goPort->go();
    //services->releasePort("goPort");
    //services->releasePort("cca.BuilderService");
  }
  catch (const sci::cca::CCAException::pointer &e) {
    std::cerr << e->getNote() << std::endl;
  }
  statusBar->SetStatusText("Components built", 0);
}

void BuilderWindow::InstantiateComponent(const sci::cca::ComponentClassDescription::pointer& cd)
{
std::cerr << "BuilderWindow::InstantiateComponent(..)" << std::endl;
  statusBar->SetStatusText("Build component", 0);
  try {
    TypeMap *tm = new TypeMap;
    tm->putString("LOADER NAME", cd->getLoaderName());

    sci::cca::ComponentID::pointer cid = builder->createInstance(cd->getComponentClassName(), sci::cca::TypeMap::pointer(tm));
    if (! cid.isNull()) {
      std::cerr << "wx: Got " << cid->getInstanceName() << std::endl;
      networkCanvas->AddIcon(cid);
    }
  }
  catch (const sci::cca::CCAException::pointer &e) {
    std::cerr << e->getNote() << std::endl;
  }
  statusBar->SetStatusText("Component built", 0);
}


///////////////////////////////////////////////////////////////////////////
// private member functions

void BuilderWindow::buildPackageMenus()
{
  //statusBar()->message("Building component menus...", 4000);
  //setCursor(Qt::WaitCursor);
  //componentMenu->clear(); // canvas popup menu
  menus.clear();

  SSIDL::array1<sci::cca::ComponentClassDescription::pointer> list;
  builder->getComponentClassDescriptions(list);

  for (std::vector<sci::cca::ComponentClassDescription::pointer>::iterator iter = list.begin(); iter != list.end(); iter++) {
    // model name could be obtained somehow locally.
    // and we can assume that the remote component model is always "CCA"
    std::string model = (*iter)->getComponentModelName();
    std::string loaderName = (*iter)->getLoaderName();

    std::string name = (*iter)->getComponentClassName();

    // component class has a loader that is not in this address space?
    if (loaderName != "") {
      std::string::size_type i = name.find_first_of(".");
      name.insert(i, "@" + loaderName);
    }
    if (menus.find(model) == menus.end()) {
      menus[model] = new MenuTree(this, url);
    }
    std::vector<std::string> splitname = split_string(name, '.');
    menus[model]->add(splitname, 0, *iter, name);
  }

  for (std::map<std::string, MenuTree*>::iterator iter = menus.begin();
       iter != menus.end(); iter++) {
    iter->second->coalesce();
  }

  //componentMenu->insertItem("Components");

  for (MenuMap::iterator iter = menus.begin(); iter != menus.end(); iter++) {
    //QPopupMenu* menu = new QPopupMenu(this);
    //menu->setFont(*bFont);
    wxMenu *menu = new wxMenu(wxT(""), wxMENU_TEAROFF);
    //int menuID;
    iter->second->populateMenu(menu);
    int menuIndex = menuBar->FindMenu(wxT(iter->first));
    // must be tested after adding components at runtime
    if (menuIndex != wxNOT_FOUND) {
      wxMenu* oldMenu = menuBar->Replace(menuIndex, menu, wxT(iter->first));
      delete oldMenu;
    } else {
      menuBar->Append(menu, wxT(iter->first));
    }
  }

  //unsetCursor();
    //statusBar()->clear();
}


}

