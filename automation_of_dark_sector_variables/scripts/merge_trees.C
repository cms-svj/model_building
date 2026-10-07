#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include <vector>
#include <string>
#include <memory>
#include <iostream>

int merge_trees(
  const char* fMain   = "events.root",
  const char* fFriend = "events_friend.root",
  const char* fOut    = "merged_tree.root",
  const char* tMain   = "Delphes",
  const char* tFr     = "DelphesFriend"
) {

  // Open inputs
  TFile fm(fMain);
  if (fm.IsZombie()) {
    std::cerr << "Cannot open " << fMain << "\n";
    return 1;
  }

  TTree* tm = (TTree*)fm.Get(tMain);

  TFile ff(fFriend);
  if (ff.IsZombie()) {
    std::cerr << "Cannot open " << fFriend << "\n";
    return 1;
  }

  TTree* tf = (TTree*)ff.Get(tFr);

  if (!tm) {
    std::cerr << "[debug] Main tree not found: " << tMain << "\n";
    std::cerr << "[debug] Keys in main file:\n";
    fm.ls();
  }

  if (!tf) {
    std::cerr << "[debug] Friend tree not found: " << tFr << "\n";
    std::cerr << "[debug] Keys in friend file:\n";
    ff.ls();
  }

  if (!tm || !tf) {
    std::cerr << "Missing tree(s).\n";
    return 2;
  }

  // --- Which friend branches to add (leaf-counted Float_t arrays) ---
  const std::vector<std::string> want = {
    "FatJet_N2", "FatJet_N3",
    "Lund_X", "Lund_Y"
  };

  // Readers
  TTreeReader rMain(tm);
  TTreeReader rFr(tf);

  // Build readers for the subset that actually exists
  struct FriendItem {
    std::string inName;                 // e.g. "FatJet_ptD"
    std::string outName;                // e.g. "F_FatJet_ptD"
    std::unique_ptr<TTreeReaderArray<Float_t>> rdr;  // leaf-counted array
    std::vector<float> buffer;          // output buffer
  };
  std::vector<FriendItem> items;
  items.reserve(want.size());
  for (const auto& nm : want) {
    if (!tf->GetBranch(nm.c_str())) {
      std::cerr << "[skip] friend branch not found: " << nm << "\n";
      continue;
    }
    FriendItem it;
    it.inName  = nm;
    it.outName = "F_" + nm;     // prefix to avoid name collisions
    it.rdr     = std::make_unique<TTreeReaderArray<Float_t>>(rFr, nm.c_str());
    items.emplace_back(std::move(it));
  }
  if (items.empty()) {
    std::cerr << "No requested friend branches found; nothing to do.\n";
    return 4;
  }

  // --- Output: clone full Delphes structure, connect addresses to tm ---
  TFile fout(fOut, "RECREATE");
  TTree* tout = tm->CloneTree(0);   // clone structure only
  tout->SetName(tMain);             // keep familiar name
  // Make tout's branches read from tm's memory so Fill() writes full Delphes
  tout->CopyAddresses(tm);

  // Add new friend branches to tout (pointing at our std::vector buffers)
  for (auto& it : items) {
    tout->Branch(it.outName.c_str(), &it.buffer);
  }

  // --- Event loop: sync main and friend, copy arrays, and fill once per event
  const Long64_t n = tm->GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    // advance both readers/trees in lockstep
    if (!rMain.Next()) break;
    if (!rFr.Next())   break;

    // copy each friend array into the vector buffer
    for (auto& it : items) {
      auto& rdr = *it.rdr;
      it.buffer.assign(rdr.begin(), rdr.end());
    }

    // Also load the main entry so tm's branch buffers are set
    tm->GetEntry(i);

    // write full Delphes + friend vectors in one go
    tout->Fill();
  }

    
  fout.cd();
  tout->Write(tMain, TObject::kOverwrite);
  fout.Close();


//   fout.Write();
  std::cout << "Wrote " << fOut << " with full Delphes + "
            << items.size() << " friend branches (prefixed F_*)\n";
  return 0;
}
