//  
//  For more information, please see: http://software.sci.utah.edu
//  
//  The MIT License
//  
//  Copyright (c) 2006 Scientific Computing and Imaging Institute,
//  University of Utah.
//  
//  License for the specific language governing rights and limitations under
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
//  
//    File   : XMLIO.cc
//    Author : McKay Davis
//    Date   : Tue Jun 27 13:01:28 2006

#include <Core/Skinner/XMLIO.h>
#include <Core/Skinner/Variables.h>
#include <Core/Skinner/Box.h>
#include <Core/Skinner/Collection.h>
#include <Core/Skinner/Grid.h>
#include <Core/Skinner/Layout.h>
#include <Core/Skinner/Text.h>
#include <Core/Skinner/Texture.h>
#include <Core/Skinner/Window.h>

#include <Core/XMLUtil/XMLUtil.h>
#include <Core/Containers/StringUtil.h>
#include <Core/Util/Environment.h>
#include <Core/Util/Assert.h>

#include <libxml/xmlreader.h>
#include <libxml/catalog.h>
#include <iostream>


namespace SCIRun {
  namespace Skinner {
    XMLIO::DrawableMakerMap_t XMLIO::makers_;

    XMLIO::XMLIO()
    {
    }

    XMLIO::~XMLIO()
    {
    }

    Drawables_t
    XMLIO::load(const string &filename)
    {
      Drawables_t objs;
      /*
       * this initialize the library and check potential ABI mismatches
       * between the version it was compiled for and the actual shared
       * library used.
       */
      
      LIBXML_TEST_VERSION;
      
      xmlParserCtxtPtr ctxt; /* the parser context */
      xmlDocPtr doc; /* the resulting document tree */
      
      string dtd = string(sci_getenv("SCIRUN_SRCDIR")) + 
        string("/Core/Skinner/skinner.dtd");
      xmlInitializeCatalog();
      xmlCatalogAdd(XMLUtil::char_to_xmlChar("public"), 
                    XMLUtil::char_to_xmlChar("-//Skinner/Drawable DTD"), 
                    XMLUtil::char_to_xmlChar(dtd.c_str()));
     
      /* create a parser context */
      ctxt = xmlNewParserCtxt();
      if (!ctxt) {
        std::cerr << "XMLIO::load failed xmlNewParserCtx()\n";
        return objs;
      }

      /* parse the file, activating the DTD validation option */
      doc = xmlCtxtReadFile(ctxt, filename.c_str(), 0, (XML_PARSE_DTDATTR | 
                                                        XML_PARSE_DTDVALID | 
                                                        XML_PARSE_PEDANTIC));
      /* check if parsing suceeded */
      if (!doc) {
        std::cerr << "Skinner::XMLIO::load failed to parse " 
                  << filename << std::endl;
        return objs;
      } 
      if (!ctxt->valid) {
          std::cerr << "Skinner::XMLIO::load dailed to validate " 
                    << filename << std::endl;
          return objs;
      }
      

      // parse the doc at network node.

      for (xmlNode *cnode=doc->children; cnode!=0; cnode=cnode->next) {
        if (XMLUtil::node_is_dtd(cnode, "skinner")) 
          continue;
        if (XMLUtil::node_is_element(cnode, "skinner")) 
          objs = eval_skinner_node(cnode, filename);
        else if (!XMLUtil::node_is_comment(cnode))
          throw "Unknown node type";
      }               

      //Drawable *obj = 0;
      xmlFreeDoc(doc);
      /* free up the parser context */
      xmlFreeParserCtxt(ctxt);  
      xmlCleanupParser();

      return objs;
    }
    


    Drawable *
    XMLIO::eval_object_node(const xmlNodePtr node, 
                            Variables *variables,
                            string_node_map_t &definitions,
                            //TargetSignalMap_t &signals,
                            SignalThrower::SignalCatchers_t &catchers) 
    {
      ASSERT(XMLUtil::node_is_element(node, "object"));

      // classname is exact class type for this skiner drawable
      string classname = XMLUtil::node_att_as_string(node, "class");

      // The classname could refer to a previously parsed <definition> node
      // in that case, we instance the contained single object node as it were 
      // createed in the current context.  <var>s and <signals> are tags
      // are allowed in both this <object> tag and the <definiition>
      // encapsulated <object> tag.  The Variables and signals of the
      // encapsulated tag are merged last.
      string_node_map_t::iterator definition = definitions.find(classname);

      // Object id is not required, create unique one if not found
      // The first merged node contains the id
      string unique_id = "";
      if (!XMLUtil::maybe_get_att_as_string(node, "id",  unique_id)) {
        // Note: This isnt absolutely guaranteed to make unique ids
        // It can be broken by hardcoding the id in the XML
        static int objcount = 0;
        unique_id = classname+"."+to_string(objcount++);
      }
      
      typedef vector<xmlNodePtr> Nodes_t;
      Nodes_t merged_nodes(1,node);

      // Is this object tag, just an alias to another object node?
      bool isa_definition = (definition != definitions.end());
      if (isa_definition) {
        // If we are just a reference to a definition, we need to switch
        // classnames to the encapsulated object
        const xmlNodePtr dnode = definition->second;
        classname = XMLUtil::node_att_as_string(dnode, "class");
        
        // When searching for vars, signals, and children
        // we need to merge the two nodes
        merged_nodes.push_back(definition->second);
      }


      // First, before we can even construct the object, we need to 
      // get the Variables that determine this instances unique properties
      // Create a new Variables context with our Unique ID
      // This creates new memory that needs to be freed by the object
      variables = variables->spawn(unique_id);

      // Spin through the xml var nodes and set their values
      for (Nodes_t::iterator mnode = merged_nodes.begin(); 
           mnode != merged_nodes.end(); ++mnode) {        
        for (xmlNode *cnode = (*mnode)->children; cnode; cnode = cnode->next) {
          if (XMLUtil::node_is_element(cnode, "var")) {
            eval_var_node(cnode, variables);
          }
        }
      }

      // Now we have Variables, Create the Object!
      Drawable * object = 0;

      // First, see if the current catchers can create and throw back 
      // an object of type "classname"
      string makerstr = classname+"_Maker";
      // The special Signal to ask for a maker is created
      event_handle_t find_maker = new MakerSignal(makerstr, variables);
      // And thrown to the Catcher...
      event_handle_t catcher_return = 
        SignalThrower::throw_signal(catchers, find_maker);
      // And we see what the cather returned
      MakerSignal *made = dynamic_cast<MakerSignal*>(catcher_return.get_rep());
      if (made && made->get_signal_name() == classname+"_Made") {
        // It returned a Maker that we wanted... Hooray!
        object = dynamic_cast<Drawable *>(made->get_signal_thrower());
      }


      // Search the static_makers table and see if it contains
      // the classname we want....DEPRECIATED, maybe goes away?
      if (makers_.find(classname) != makers_.end()) {
        object = (*makers_[classname])(variables);
      }
      

      // At this point, the if the object is uninstantiatable, return
      if (!object) {
        delete variables;
        cerr << "Skinner::XMLIO::eval_object_node - object class=\"" 
             << classname << "\" cannot find maker\n";
        return 0;
      }


      // Now the object is created, fill the catchersondeck tree
      // with targets contained by the new object.
      cerr << object->get_id() << " - adding catcher";
      SignalCatcher::TargetIDs_t catcher_targets =object->get_all_target_ids();
      SignalCatcher::TargetIDs_t::iterator id_iter;
      for(id_iter  = catcher_targets.begin();
          id_iter != catcher_targets.end(); ++id_iter) {
        const string &id = *id_iter;
        SignalThrower::Catcher_t catcher_target = 
          make_pair(object, object->get_catcher(id));
        catchers[id].push_front(catcher_target);
      }


      // Now the Catchers On Deck are ready, look if the xml file has
      // created any signals to hookup
      for (Nodes_t::iterator mnode = merged_nodes.begin();
           mnode != merged_nodes.end(); ++mnode) {
        for (xmlNode *cnode = (*mnode)->children; cnode; cnode = cnode->next) {
          if (XMLUtil::node_is_element(cnode, "signal")) {
            eval_signal_node(cnode, object, catchers);
          } 
        }
      }
      
      // Time to look for children object nodes
      Drawables_t children(0);
      bool have_unwanted_children = false;
      Parent *parent = dynamic_cast<Parent *>(object);
      // Search the merged nodes for object nodes.
      for (Nodes_t::iterator mnode = merged_nodes.begin(); 
           mnode != merged_nodes.end(); ++mnode) {        
        for (xmlNode *cnode = (*mnode)->children; cnode; cnode = cnode->next) {
          if (XMLUtil::node_is_element(cnode, "object")) {
            if (parent) { 
              Drawable *child = 
                eval_object_node(cnode, variables, definitions, catchers);
              if (child) {
                children.push_back(child);
              }
            } else {
              have_unwanted_children = true;
            }
          }
        }
      }
      
      if (parent && children.size()) {
        parent->set_children(children);
      }
            
      if (have_unwanted_children) { 
          cerr << "class : " << classname << " does not allow <object>\n";
      }

      // Re-get all object ids, as we may have pushed some aliases 
      // during eval_skinner_node
      cerr << object->get_id() << " - removing catcher";
      catcher_targets = object->get_all_target_ids();
      for(id_iter  = catcher_targets.begin();
          id_iter != catcher_targets.end(); ++id_iter) {
        const string &id = *id_iter;
        ASSERTMSG(!catchers[id].empty(), "Catchers Empty!");
        catchers[id].pop_front();
      }
      
      return object;
    }


    vector<Drawable *>
    XMLIO::eval_skinner_node(const xmlNodePtr node, const string &id)
    {
      ASSERT(XMLUtil::node_is_element(node, "skinner"));
      Drawables_t children;
      string_node_map_t definitions;
      //      TargetSignalMap_t signals;
      SignalThrower::SignalCatchers_t catchers;
      Variables *variables = new Variables(id);

      for (xmlNode *cnode=node->children; cnode!=0; cnode=cnode->next) {
        if (XMLUtil::node_is_element(cnode, "definition")) {
          eval_definition_node(cnode, definitions);
        } else if (XMLUtil::node_is_element(cnode, "object")) {
          children.push_back(eval_object_node(cnode, 
                                              variables, 
                                              definitions,
                                              catchers));
        } 
      }
      ASSERT(!children.empty());
      return children;
    }

    void
    XMLIO::eval_definition_node(const xmlNodePtr node,
                                string_node_map_t &definitions)
    {
      ASSERT(XMLUtil::node_is_element(node, "definition"));
      string classname = XMLUtil::node_att_as_string(node, "class");

      for (xmlNode *cnode=node->children; cnode!=0; cnode=cnode->next) {
        if (XMLUtil::node_is_element(cnode, "object")) {
          definitions[classname] = cnode;
        } 
      }
    }

    void
    XMLIO::eval_var_node(const xmlNodePtr node,
                         Variables *variables) 
    {
      ASSERT(XMLUtil::node_is_element(node, "var"));
      ASSERT(variables);
      bool propagate = 
        XMLUtil::node_att_as_string(node, "propagate") == "yes" ? true : false;

      variables->insert(XMLUtil::node_att_as_string(node, "name"),
                        XMLUtil::xmlChar_to_char(node->children->content),
                        XMLUtil::node_att_as_string(node, "type"),
                        propagate);
    }
                        



    void
    XMLIO::eval_signal_node(const xmlNodePtr node,
                            Drawable *object,
                            SignalThrower::SignalCatchers_t &catchers)
      //                            TargetSignalMap_t &signals) 
    {
      ASSERT(XMLUtil::node_is_element(node, "signal"));
      string signalname = XMLUtil::node_att_as_string(node, "name");
      string signaltarget = XMLUtil::node_att_as_string(node, "target");
      SignalThrower *thrower = dynamic_cast<SignalThrower *>(object);
      if (!thrower) {
        cerr << "Signal in non throwing object\n";
        return;
      }


      SignalThrower::SignalCatchers_t::iterator catchers_iter = 
        catchers.find(signaltarget);

      if (catchers_iter == catchers.end()) {
        cerr << "Signal " << signalname 
             << " cannot find target " << signaltarget << std::endl;
        return;
      }

      SignalCatcher *thrower_as_catcher = dynamic_cast<SignalCatcher*>(object);


//        cerr << "Catchers for " << object->get_id() << ": ";
//       SignalThrower::SignalCatchers_t::iterator citer = catchers.begin();
//       for (; citer != catchers.end(); ++citer) {
//         cerr << citer->first << (citer->second.empty() ? 0 : 1) << ", ";
//       }
//       cerr << std::endl;

//       if(catchers.find(signaltarget) == catchers.end()) {
//         cerr << "Signal " << signalname 
//              << " cannot find target " << signaltarget << std::endl;
//         return;
//       }

      SignalThrower::CatchersOnDeck_t &signal_catchers = catchers_iter->second;
      SignalThrower::CatchersOnDeck_t::iterator catcher_iter = 
        signal_catchers.begin();
      for (;catcher_iter != signal_catchers.end(); ++catcher_iter) {
        SignalCatcher *catcher = catcher_iter->first;
        SignalCatcher::CatcherFunctionPtr function = catcher_iter->second;
        //          catcher->get_catcher(signaltarget);
        if (function) {
          if (thrower->get_signal_id(signalname)) {
            cerr << " Add Catcher: " << signalname << std::endl;
            thrower->catchers_[signalname].push_front
              (make_pair(catcher, function));
          }  else {
            cerr << "Translating signal: " 
                 << signalname << " to: " << signaltarget << std::endl;
            catchers[signalname].push_front(make_pair(catcher, function));
            thrower_as_catcher->catcher_functions_[signalname] = function;
            cerr << object->get_id() << " - now contains";
            object->get_all_target_ids();
          }


//           if (thrower->add_catcher_function(signalname, catcher, function)) {
//             Drawable *catcherdrawable = static_cast<Drawable *>(catcher);
//             if (catcherdrawable) {
//               cerr << "Catcher for " << signalname 
//                    << " is id " << catcherdrawable->get_id() 
//                    << " target " << signaltarget << std::endl;
//             }
//           } else if (thrower_as_catcher) {            
//             cerr << "Translating signal: " 
//                  << signalname << " to: " << signaltarget << std::endl;
//             catchers[signalname].push_front(make_pair(catcher, function));
//             thrower_as_catcher->catcher_functions_[signalname] = function;
//           }
        }
      }
    }
        
        
      


      //      for (SignalCatcherList_t::iterator catcher_iter = catchers.begin();
//              catcher_iter != catchers.end(); ++catcher_iter) {
//           SignalCatcher *catcher = *catcher_iter;
//           SignalCatcher::CatcherFunctionPtr function = 
//             catcher->get_catcher(signaltarget);
//           if (function) {
//             thrower->add_catcher_function(signalname, catcher, function);
//           }
//         }
//       }
      


//        int signal_id = thrower->get_signal_id(signalname);
//       if (signal_id) {
//         signals[signaltarget] = make_pair(signalname, thrower);
//       }


//         for (SignalCatcherList_t::iterator catcher_iter = catchers.begin();
//              catcher_iter != catchers.end(); ++catcher_iter) {
//           SignalCatcher *catcher = *catcher_iter;
//           SignalCatcher::CatcherFunctionPtr function = 
//             catcher->get_catcher(signaltarget);
//           if (function) {
//             thrower->add_catcher_function(signalname, catcher, function);
//           }
//         }
//       }
//    }
                        
      
      
    

    void
    XMLIO::register_maker(const string &name, DrawableMakerFunc_t *maker)
    {
      makers_[name] = maker;
    }





  }
}
