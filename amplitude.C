// Mainframe macro generated from application: /home/cami20/root/bin/root.exe
// By ROOT version 6.28/06 on 2025-05-20 09:44:49

#ifndef ROOT_TGFrame
#include "TGFrame.h"
#endif
#ifndef ROOT_TGListBox
#include "TGListBox.h"
#endif
#ifndef ROOT_TGScrollBar
#include "TGScrollBar.h"
#endif
#ifndef ROOT_TGComboBox
#include "TGComboBox.h"
#endif
#ifndef ROOT_TGMenu
#include "TGMenu.h"
#endif
#ifndef ROOT_TGFileDialog
#include "TGFileDialog.h"
#endif
#ifndef ROOT_TGButtonGroup
#include "TGButtonGroup.h"
#endif
#ifndef ROOT_TGCanvas
#include "TGCanvas.h"
#endif
#ifndef ROOT_TGFSContainer
#include "TGFSContainer.h"
#endif
#ifndef ROOT_TGButton
#include "TGButton.h"
#endif
#ifndef ROOT_TRootContextMenu
#include "TRootContextMenu.h"
#endif
#ifndef ROOT_TGFSComboBox
#include "TGFSComboBox.h"
#endif
#ifndef ROOT_TGLabel
#include "TGLabel.h"
#endif
#ifndef ROOT_TGListView
#include "TGListView.h"
#endif
#ifndef ROOT_TGSplitter
#include "TGSplitter.h"
#endif
#ifndef ROOT_TGTextEntry
#include "TGTextEntry.h"
#endif
#ifndef ROOT_TRootCanvas
#include "TRootCanvas.h"
#endif
#ifndef ROOT_TGDockableFrame
#include "TGDockableFrame.h"
#endif
#ifndef ROOT_TG3DLine
#include "TG3DLine.h"
#endif
#ifndef ROOT_TGToolTip
#include "TGToolTip.h"
#endif

#include "Riostream.h"

void amplitude()
{

   // main frame
   TGMainFrame *fRootCanvas2 = new TGMainFrame(gClient->GetRoot(),10,10,kMainFrame | kVerticalFrame);

   // menu bar
   TGMenuBar *fMenuBar12 = new TGMenuBar(fRootCanvas2,1908,22,kHorizontalFrame);

   // "File" menu
   TGPopupMenu *fPopupMenu4 = new TGPopupMenu(gClient->GetDefaultRoot(),119,144,kOwnBackground);
   fPopupMenu4->AddEntry("&New Canvas",0);
   fPopupMenu4->AddEntry("&Open...",1);
   fPopupMenu4->AddEntry("&Close Canvas",13);
   fPopupMenu4->AddSeparator();

   // cascaded menu "Save"
   TGPopupMenu *fPopupMenu3 = new TGPopupMenu(gClient->GetDefaultRoot(),70,168,kOwnBackground);
   fPopupMenu3->AddEntry("c.&ps",5);
   fPopupMenu3->AddEntry("c.&eps",6);
   fPopupMenu3->AddEntry("c.p&df",7);
   fPopupMenu3->AddEntry("c.&tex",11);
   fPopupMenu3->AddEntry("c.&gif",8);
   fPopupMenu3->AddEntry("c.&jpg",9);
   fPopupMenu3->AddEntry("c.&png",10);
   fPopupMenu3->AddEntry("c.&C",4);
   fPopupMenu3->AddEntry("c.&root",3);
   fPopupMenu4->AddPopup("&Save",fPopupMenu3);
   fPopupMenu4->AddEntry("Save &As...",2);
   fPopupMenu4->AddSeparator();
   fPopupMenu4->AddEntry("&Print...",12);
   fPopupMenu4->AddSeparator();
   fPopupMenu4->AddEntry("&Quit ROOT",14);
   fMenuBar12->AddPopup("&File",fPopupMenu4, new TGLayoutHints(kLHintsLeft | kLHintsTop,0,4,0,0));

   // "Edit" menu
   TGPopupMenu *fPopupMenu6 = new TGPopupMenu(gClient->GetDefaultRoot(),80,144,kOwnBackground);
   fPopupMenu6->AddEntry("&Style...",15);
   fPopupMenu6->AddSeparator();
   fPopupMenu6->AddEntry("Cu&t",16);
   fPopupMenu6->DisableEntry(16);
   fPopupMenu6->AddEntry("&Copy",17);
   fPopupMenu6->DisableEntry(17);
   fPopupMenu6->AddEntry("&Paste",18);
   fPopupMenu6->DisableEntry(18);
   fPopupMenu6->AddSeparator();

   // cascaded menu "Clear"
   TGPopupMenu *fPopupMenu5 = new TGPopupMenu(gClient->GetDefaultRoot(),79,42,kOwnBackground);
   fPopupMenu5->AddEntry("&Pad",19);
   fPopupMenu5->AddEntry("&Canvas",20);
   fPopupMenu6->AddPopup("C&lear",fPopupMenu5);
   fPopupMenu6->AddSeparator();
   fPopupMenu6->AddEntry("&Undo",21);
   fPopupMenu6->DisableEntry(21);
   fPopupMenu6->AddEntry("&Redo",22);
   fPopupMenu6->DisableEntry(22);
   fMenuBar12->AddPopup("&Edit",fPopupMenu6, new TGLayoutHints(kLHintsLeft | kLHintsTop,0,4,0,0));

   // "View" menu
   TGPopupMenu *fPopupMenu8 = new TGPopupMenu(gClient->GetDefaultRoot(),141,180,kOwnBackground);
   fPopupMenu8->AddEntry("&Editor",23);
   fPopupMenu8->AddEntry("&Toolbar",24);
   fPopupMenu8->AddEntry("Event &Statusbar",25);
   fPopupMenu8->AddEntry("T&oolTip Info",26);
   fPopupMenu8->AddSeparator();
   fPopupMenu8->AddEntry("&Colors",27);
   fPopupMenu8->AddEntry("&Fonts",28);
   fPopupMenu8->DisableEntry(28);
   fPopupMenu8->AddEntry("&Markers",29);
   fPopupMenu8->AddSeparator();
   fPopupMenu8->AddEntry("&Iconify",30);
   fPopupMenu8->AddSeparator();

   // cascaded menu "View With"
   TGPopupMenu *fPopupMenu7 = new TGPopupMenu(gClient->GetDefaultRoot(),84,42,kOwnBackground);
   fPopupMenu7->AddEntry("&X3D",31);
   fPopupMenu7->AddEntry("&OpenGL",32);
   fPopupMenu8->AddPopup("&View With",fPopupMenu7);
   fMenuBar12->AddPopup("&View",fPopupMenu8, new TGLayoutHints(kLHintsLeft | kLHintsTop,0,4,0,0));

   // "Options" menu
   TGPopupMenu *fPopupMenu9 = new TGPopupMenu(gClient->GetDefaultRoot(),169,216,kOwnBackground);
   fPopupMenu9->AddEntry("&Auto Resize Canvas",33);
   fPopupMenu9->CheckEntry(33);
   fPopupMenu9->AddEntry("&Resize Canvas",34);
   fPopupMenu9->AddEntry("&Move Opaque",35);
   fPopupMenu9->CheckEntry(35);
   fPopupMenu9->AddEntry("Resize &Opaque",36);
   fPopupMenu9->CheckEntry(36);
   fPopupMenu9->AddSeparator();
   fPopupMenu9->AddEntry("&Interrupt",37);
   fPopupMenu9->AddEntry("R&efresh",38);
   fPopupMenu9->AddSeparator();
   fPopupMenu9->AddEntry("&Pad Auto Exec",39);
   fPopupMenu9->AddSeparator();
   fPopupMenu9->AddEntry("&Statistics",40);
   fPopupMenu9->CheckEntry(40);
   fPopupMenu9->AddEntry("Histogram &Title",41);
   fPopupMenu9->CheckEntry(41);
   fPopupMenu9->AddEntry("&Fit Parameters",42);
   fPopupMenu9->AddEntry("Can Edit &Histograms",43);
   fMenuBar12->AddPopup("&Options",fPopupMenu9, new TGLayoutHints(kLHintsLeft | kLHintsTop,0,4,0,0));

   // "Tools" menu
   TGPopupMenu *fPopupMenu10 = new TGPopupMenu(gClient->GetDefaultRoot(),134,114,kOwnBackground);
   fPopupMenu10->AddEntry("&Inspect ROOT",44);
   fPopupMenu10->AddEntry("&Class Tree",45);
   fPopupMenu10->AddEntry("&Fit Panel",46);
   fPopupMenu10->AddEntry("&Start Browser",47);
   fPopupMenu10->AddEntry("&Gui Builder",48);
   fPopupMenu10->AddEntry("&Event Recorder",49);
   fMenuBar12->AddPopup("&Tools",fPopupMenu10, new TGLayoutHints(kLHintsLeft | kLHintsTop,0,4,0,0));

   // "Help" menu
   TGPopupMenu *fPopupMenu11 = new TGPopupMenu(gClient->GetDefaultRoot(),138,158,kOwnBackground);
   fPopupMenu11->AddLabel("Basic Help On...");
   fPopupMenu11->DefaultEntry(-1);
   fPopupMenu11->AddSeparator();
   fPopupMenu11->AddEntry("&Canvas",51);
   fPopupMenu11->AddEntry("&Menus",52);
   fPopupMenu11->AddEntry("&Graphics Editor",53);
   fPopupMenu11->AddEntry("&Browser",54);
   fPopupMenu11->AddEntry("&Objects",55);
   fPopupMenu11->AddEntry("&PostScript",56);
   fPopupMenu11->AddSeparator();
   fPopupMenu11->AddEntry("&About ROOT...",50);
   fMenuBar12->AddPopup("&Help",fPopupMenu11, new TGLayoutHints(kLHintsRight | kLHintsTop));
   fRootCanvas2->AddFrame(fMenuBar12, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,0,0,1,1));
   TGHorizontal3DLine *fHorizontal3DLine20 = new TGHorizontal3DLine(fRootCanvas2,800,2);
   fRootCanvas2->AddFrame(fHorizontal3DLine20, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

   // dockable frame
   TGDockableFrame *fDockableFrame21 = new TGDockableFrame(fRootCanvas2);

   // next lines belong to the dockable frame widget
   fDockableFrame21->EnableUndock(kTRUE);
   fDockableFrame21->EnableHide(kFALSE);
   fDockableFrame21->SetWindowName("ToolBar: Crossover Analysis");
   fDockableFrame21->DockContainer();

   fRootCanvas2->AddFrame(fDockableFrame21, new TGLayoutHints(kLHintsExpandX));
   TGHorizontal3DLine *fHorizontal3DLine26 = new TGHorizontal3DLine(fRootCanvas2,800,2);
   fRootCanvas2->AddFrame(fHorizontal3DLine26, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

   // composite frame
   TGCompositeFrame *fCompositeFrame27 = new TGCompositeFrame(fRootCanvas2,1908,963,kHorizontalFrame);

   // composite frame
   TGCompositeFrame *fCompositeFrame28 = new TGCompositeFrame(fCompositeFrame27,175,576,kFixedWidth);
   fCompositeFrame28->SetLayoutManager(new TGVerticalLayout(fCompositeFrame28));

   fCompositeFrame27->AddFrame(fCompositeFrame28, new TGLayoutHints(kLHintsLeft | kLHintsExpandY));

   // canvas widget
   TGCanvas *fCanvas29 = new TGCanvas(fCompositeFrame27,1908,963,kSunkenFrame);

   // canvas viewport
   TGViewPort *fViewPort30 = fCanvas29->GetViewPort();

   // canvas container
   Int_t canvasID = gVirtualX->InitWindow((ULongptr_t)fCanvas29->GetId());
   Window_t winC = gVirtualX->GetWindowID(canvasID);
   TGCompositeFrame *fCompositeFrame39 = new TGCompositeFrame(gClient,winC,fViewPort30);
   fViewPort30->AddFrame(fCompositeFrame39);
   fCompositeFrame39->SetLayoutManager(new TGVerticalLayout(fCompositeFrame39));
   fCompositeFrame39->MapSubwindows();
   fCanvas29->SetContainer(fCompositeFrame39);
   fCanvas29->MapSubwindows();
   fCompositeFrame27->AddFrame(fCanvas29, new TGLayoutHints(kLHintsRight | kLHintsExpandX | kLHintsExpandY));

   fRootCanvas2->AddFrame(fCompositeFrame27, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   // status bar
   TGStatusBar *fStatusBar42 = new TGStatusBar(fRootCanvas2,10,18);
   Int_t partsusBar42[] = {33,10,10,47};
   fStatusBar42->SetParts(partsusBar42,4);
   fRootCanvas2->AddFrame(fStatusBar42, new TGLayoutHints(kLHintsLeft | kLHintsBottom | kLHintsExpandX,2,2,1,1));

   fRootCanvas2->SetWindowName("Crossover Analysis");
   fRootCanvas2->SetIconName("Crossover Analysis");
   fRootCanvas2->SetIconPixmap("macro_s.xpm");
   fRootCanvas2->SetClassHints("ROOT","Canvas");
   fRootCanvas2->SetMWMHints(kMWMDecorAll,
                        kMWMFuncAll,
                        kMWMInputModeless);
   fRootCanvas2->MapSubwindows();
   fHorizontal3DLine20->UnmapWindow();
   fDockableFrame21->UnmapWindow();
   fHorizontal3DLine26->UnmapWindow();
   fCompositeFrame28->UnmapWindow();
   fStatusBar42->UnmapWindow();

   fRootCanvas2->Resize(fRootCanvas2->GetDefaultSize());
   fRootCanvas2->MapWindow();
   fRootCanvas2->Resize(1908,987);
}  
