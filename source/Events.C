#define Events_cxx
#include "Events.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>
#include<cmath>
#include "TLegend.h"
#include "TFractionFitter.h"

int nev=0; // number of events
int nmu=0; // number of muons
int nele=0; //number of electrons
int ntau=0; // number of taus



//FUNCTION THAT RETURNS THE INDEX OF THE PARTICLE; IF THAT PARTICLE IS A LEPTON

Int_t Events::TypeIndex()
{
    for (int i=0; i<=nLHEPart; i++) { //loop over the particles
        
        if(abs(LHEPart_pdgId[i])==11||abs(LHEPart_pdgId[i])==13||abs(LHEPart_pdgId[i])==15)  return i; //if a particle is a lepton, return the index
    }
    return -1; //else return -1
}

//FUNCTION THAT FILLS THE HISTOGRAMS GIVEN

void Events::Filling(Int_t ind, Int_t ilept)
{
    
    ((h[ind]))->Fill(LHEPart_pt[ilept]); //fill histogram for pt WHITOUT CUTS
    
    if(nElectron>0){ //if there is at least one reconstructed Electron
        ((hR[2*ind]))->Fill(Electron_pt[0]); //store pT value of Electron
        
        hPlog[2*ind]->Fill(std::log(std::abs(Electron_dxy[0])));
        
        if(Electron_pfRelIso03_all[0]<0.15&&Electron_cutBased[0]==4){ //cuts on the cone, i want an isolated particle
            ((hRi[2*ind]))->Fill(Electron_pt[0]); //fill histogram
            
            
            hP[2*ind]->Fill(std::abs(Electron_dxy[0]));
            
            hPz[2*ind]->Fill(std::abs(Electron_dz[0]));
            
            if(Electron_dxyErr[0]!=0)
            hPerr[2*ind]->Fill(std::abs(Electron_dxy[0])/Electron_dxyErr[0]);

            if(Electron_dzErr[0]!=0)
            hPzerr[2*ind]->Fill(std::abs(Electron_dz[0])/Electron_dzErr[0]);
            
            //RctIsoPt[2*ind]=Electron_pt[0]; //store value for tree
            //top->Fill(); //fill tree
        }
    }
    if (nMuon>0){ //Muon, structure similar as above
        ((hR[2*ind+1]))->Fill(Muon_pt[0]);
        
        
        if(Muon_pfRelIso03_all[0]<0.15&&Muon_tightId[0]==true){
            ((hRi[2*ind+1]))->Fill(Muon_pt[0]);
            //RctIsoPt[2*ind+1]=Muon_pt[0];
            hPlog[2*ind+1]->Fill(std::log(std::abs(Muon_dxy[0])));
            
            hP[2*ind+1]->Fill(std::abs(Muon_dxy[0]));

            if(ind>0)data->Fill(std::abs(Muon_dxy[0]));
            
            
            hPz[2*ind+1]->Fill(std::abs(Muon_dz[0]));
            
            if(Muon_dxyErr[0]!=0)
            hPerr[2*ind+1]->Fill(std::abs(Muon_dxy[0])/Muon_dxyErr[0]);

            if(Muon_dzErr[0]!=0)
            hPzerr[2*ind+1]->Fill(std::abs(Muon_dz[0])/Muon_dzErr[0]);
            
            //top->Fill();
        }
    }
    
    return;
}

//generate three plots, one for lept_type (event type), containing the Reconstructed (and isolated) Electron and Muon pT, and the two stacked together


//e dal tau hanno meno energia

// se l'avessimo ricostruti ci sono molto più eletrroni a impulso trasverso che mancano
//c'è un problema di ricotruzione che manca. Una grossa fetta di questi elettroni non riescono ad essere ricostruiti. Me ne aspetterei 6000 ma ne ho 3000. Mi manca quello a bassissimo pT. Se aggiunto mi porterebbe la media 20. Parto da 60 e lo distribuisco su tre corpi (un eletrrone e due neutrini).

//I muoni un po' meglio. Si ricostruiscono meeglio. Hanno lo spetro più soffice (55 anziché 61) perché si ricostruiscono meglio quelli a basso pt. quindi ancora meglio per il tau. Infatti si vede meglio la zona sullo zero.

//i tau vanno in eletroni col 18%

//Nelle legende mettiamoci i numeri, pT medio

//la parte verde il 10% in più ricostruisce solo.

//la normalizzazione viene fatta solo ad alto pT. Quando vedo misure sotto 50 vedo il contributo del tau.

//Metto le slide

//Ci sono diversi tagli. Qualche indicazione c'è sull'articolo. Sono collegati all'identificazione. Ad esempio l'elettrone ha traccia e lascia la segnatura nel calorimetro elettromagnetico. Lo fanno solo i fotoni e gli elettroni. Avendo visto la traccia escludo i fotoni. E/p energia nel calorimentro e impulso misurato nel traccaitore dalla curvatura del campo magnetico. Se è un elttrone voglio E/p=1.  23 lunghezze di radiazione corrispondono una lunghezza di interazione nucleare (pione lascia poca energia). lo sciame adronico è limitato.
//Ci sono altre varibili che corrispondono allo sciame. Sono contentui dentro il raggio di moliere. Se arriva un pione lo sciame è molto più largo.
//
//con i muoni si può usare anche medium.


/*

 

*/
//Muoni, ho la traccia, non ho energia nel calorimetro. Faccio corrispondere gli hit nelle camere a muoni con la traccia.

//Guarderemo i parametri d'impatto (quanto la traccia si avvicina al vertice primario). Oltre ad avere un pt diverso, i diretti dal W puntano all'interazione principale (circa 0). La componente dal tau avrà un parametro d'impatto più grande.

//generate one additional plot, pT of particle (divided by event type) without cuts stacked together

void Events::ReconStack()

//COLOR THE HISTO GIVEN THE PARTICLE
{
    //COLOR THE HISTOGRAM
    for (int i=0; i<5;  i+=2) { //Electron Blue
        ((hR[i]))->SetFillColor(38);
        ((hRi[i]))->SetFillColor(38);
    }
    ((h[0]))->SetFillColor(38); //without cuts
    
    for (int i=1; i<6;  i+=2) { //Muon Red
        ((hR[i]))->SetFillColor(45);
        ((hRi[i]))->SetFillColor(45);
    }
    (h[1])->SetFillColor(45);
    
    (h[2])->SetFillColor(30); //Tau Green


    
    
 //STACK ALL THE HISTO

    //STACK ALL THE HISTOGRAMS
    for (int m=0; m<3; m++) { //PT STACKING IN HS[0]
        //h[m]->SetStats(kFALSE);

        (hs[0])->Add(h[m]);
    }
    
    for (int i=1,j=0; i<4; i++,j+=2) { // RECONSTROCTUD STACKING HS[1,2,3]
        ((hs[i]))->Add((hR[j]));
        ((hs[i]))->Add((hR[j+1]));
        ((hs[i+3]))->Add((hRi[j])); //AND ISO REC HS[4,5,6]
        ((hs[i+3]))->Add((hRi[j+1]));
    }
    
// DRAW THE RECON HISTOGRAM
    
  



    //DRAW JUST THE PT STACKED, WITHOUT ANY CUTS
        (c[0])->Divide(2,2);
    
    for (int i=0; i<3; i++) {
        c[0]->cd(i+1);
        h[i]->Draw();
        h[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
        h[i]->GetYaxis()->SetTitle("Events");
        
    }
        
        (c[0])->cd(4);

    
        (hs[0])->Draw();

    auto legend_wo = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend_wo->AddEntry(h[0],"Electron Events","f");
    legend_wo->AddEntry(h[1],"Muon Events","f");
    legend_wo->AddEntry(h[2],"Tau Events","f");
    legend_wo->Draw();
    
    
    //Draw stacked histograms iso and rec particles
    
    for (int i=1,j=0; i<4; i++,j+=2) {
        ((c[i]))->Divide(3,2);
        
        ((c[i]))->cd(1);
        ((hs[i]))->Draw();
        hs[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
        hs[i]->GetYaxis()->SetTitle("Events");
        auto legend_ele1 = new TLegend(0.6, 0.7, 0.9, 0.9);
        legend_ele1->AddEntry(hR[j],"Electron pT","f");
        legend_ele1->AddEntry(hR[j+1],"Muon pT","f");
        legend_ele1->Draw();
        
        ((c[i]))->cd(2);
        ((hR[j]))->Draw();
        hR[j]->GetXaxis()->SetTitle("p_{T} [GeV]");
        hR[j]->GetYaxis()->SetTitle("Events");
        
        ((c[i]))->cd(3);
        ((hR[j+1]))->Draw();
        hR[j+1]->GetXaxis()->SetTitle("p_{T} [GeV]");
        hR[j+1]->GetYaxis()->SetTitle("Events");
        
        ((c[i]))->cd(4);
        ((hs[i+3]))->Draw();
        hs[i+3]->GetXaxis()->SetTitle("p_{T} [GeV]");
        hs[i+3]->GetYaxis()->SetTitle("Events");
        auto legend_ele2 = new TLegend(0.6, 0.7, 0.9, 0.9);
        legend_ele2->AddEntry(hRi[j],"Electron pT","f");
        legend_ele2->AddEntry(hRi[j+1],"Muon pT","f");
        legend_ele2->Draw();
        
        ((c[i]))->cd(5);
        ((hRi[j]))->Draw();
        hRi[j]->GetXaxis()->SetTitle("p_{T} [GeV]");
        hRi[j]->GetYaxis()->SetTitle("Events");
        
        ((c[i]))->cd(6);
        ((hRi[j+1]))->Draw();
        hRi[j+1]->GetXaxis()->SetTitle("p_{T} [GeV]");
        hRi[j+1]->GetYaxis()->SetTitle("Events");

        c[i]->Modified();
    }
/*
    c[1]->cd(1);
    
    auto legend_ele1 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend_ele1->AddEntry(hR[0],"Electron pT","f");
    legend_ele1->AddEntry(hR[1],"Muon pT","f");
    legend_ele1->Draw();
    
    c[1]->cd(4);
    
    auto legend_ele2 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend_ele2->AddEntry(hRi[0],"Electron pT","f");
    legend_ele2->AddEntry(hRi[1],"Muon pT","f");
    legend_ele2->Draw();
    
    c[2]->cd(1);
    
    auto legend_mu1 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend_mu1->AddEntry(hR[2],"Electron pT","f");
    legend_mu1->AddEntry(hR[3],"Muon pT","f");
    legend_mu1->Draw();
    
    c[2]->cd(4);
    
    auto legend_ele2 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend_ele2->AddEntry(hRi[0],"Electron pT","f");
    legend_ele2->AddEntry(hRi[1],"Muon pT","f");
    legend_ele2->Draw();
    */
    
    c[0]->Modified();
    
     
    c[0]->SaveAs("./Graphs/new/pT_stacked_woCuts.pdf");
    c[1]->SaveAs("./Graphs/new/pT_stacked_RctIso_Ele.pdf");
    c[2]->SaveAs("./Graphs/new/pT_stacked_RctIso_Muon.pdf");
    c[3]->SaveAs("./Graphs/new/pT_stacked_RctIso_Tau.pdf");
    
    

    return;
}

//Draw two histograms. One for reconstructed isolated Electron pT, stacked by event type (lept_type), one for reconstructed muon pT, stacked by event type (lept type)

void Events::IsoStack()
{
    
    for (int i=0; i<6; i++) {
        hRi[i]->SetMarkerColor(1);
        hP[i]->SetMarkerColor(1);
        
        hPerr[i]->SetMarkerColor(1);
        
        hPz[i]->SetMarkerColor(1);
        
        hPzerr[i]->SetMarkerColor(1);
        
        hPlog[i]->SetMarkerColor(1);
    }
    
    //RECOLOR THE HISTOGRAM BY LEPT_TYPE
    ((hRi[0]))->SetFillColor(38); //electron blue
    ((hRi[1]))->SetFillColor(38);
    ((hP[0]))->SetFillColor(38);
    ((hPerr[0]))->SetFillColor(38);
    ((hP[1]))->SetFillColor(38);
    ((hPerr[1]))->SetFillColor(38);
    
    ((hPz[0]))->SetFillColor(38);
    ((hPzerr[0]))->SetFillColor(38);
    ((hPz[1]))->SetFillColor(38);
    ((hPzerr[1]))->SetFillColor(38);
    
    ((hPlog[0]))->SetFillColor(38);
    ((hPlog[1]))->SetFillColor(38);
    
    ((hRi[2]))->SetFillColor(45); //muon red
    ((hRi[3]))->SetFillColor(45);
    ((hP[2]))->SetFillColor(45);
    ((hPerr[2]))->SetFillColor(45);
    ((hP[3]))->SetFillColor(45);
    ((hPerr[3]))->SetFillColor(45);
    
    ((hPlog[2]))->SetFillColor(45);
    ((hPlog[3]))->SetFillColor(45);
    
    ((hPz[2]))->SetFillColor(45);
    ((hPzerr[2]))->SetFillColor(45);
    ((hPz[3]))->SetFillColor(45);
    ((hPzerr[3]))->SetFillColor(45);
    
    ((hRi[4]))->SetFillColor(30);  //tau red
    ((hRi[5]))->SetFillColor(30);
    ((hP[4]))->SetFillColor(30);
    ((hPerr[4]))->SetFillColor(30);
    ((hP[5]))->SetFillColor(30);
    ((hPerr[5]))->SetFillColor(30);
    
    ((hPlog[4]))->SetFillColor(30);
    ((hPlog[5]))->SetFillColor(30);
    
    ((hPz[4]))->SetFillColor(30);
    ((hPzerr[4]))->SetFillColor(30);
    ((hPz[5]))->SetFillColor(30);
    ((hPzerr[5]))->SetFillColor(30);
    
    //REFILL
    
    // RECONSTRUCTED ISOLATED STACKING OVER LEPT_TYPE, ONE FOR ELECTRON PT, THE OTHER FOR MUON PT
    for (int j=0; j<5; j+=2) {
        hRi[j]->SetStats(kFALSE); // no statistics box
        hRi[j+1]->SetStats(kFALSE); // no statistics box
        
        hP[j]->SetStats(kFALSE); // no statistics box
        hP[j+1]->SetStats(kFALSE); // no statistics box
        hPerr[j]->SetStats(kFALSE); // no statistics box
        hPerr[j+1]->SetStats(kFALSE); // no statistics box
        
        hPlog[j]->SetStats(kFALSE); // no statistics box
        hPlog[j+1]->SetStats(kFALSE); // no statistics box
        
        hPz[j]->SetStats(kFALSE); // no statistics box
        hPz[j+1]->SetStats(kFALSE); // no statistics box
        hPzerr[j]->SetStats(kFALSE); // no statistics box
        hPzerr[j+1]->SetStats(kFALSE); // no statistics box
        
        ((hs[7]))->Add((hRi[j]));
        ((hs[8]))->Add((hRi[j+1]));
        //(hs[9])->Add(hP[j]);
        //(hs[11])->Add(hPerr[j]);
        //(hs[10])->Add(hP[j+1]);
        //(hs[12])->Add(hPerr[j+1]);

    }
    
    (hs[9])->Add(hP[0]);
    (hs[9])->Add(hP[4]);
    (hs[11])->Add(hPerr[0]);
    (hs[11])->Add(hPerr[4]);
    (hs[10])->Add(hP[3]);
    (hs[10])->Add(hP[5]);
    (hs[12])->Add(hPerr[3]);
    (hs[12])->Add(hPerr[5]);
    
    (hs[13])->Add(hPz[0]);
    (hs[13])->Add(hPz[4]);
    (hs[15])->Add(hPzerr[0]);
    (hs[15])->Add(hPzerr[4]);
    (hs[14])->Add(hPz[3]);
    (hs[14])->Add(hPz[5]);
    (hs[16])->Add(hPzerr[3]);
    (hs[16])->Add(hPzerr[5]);
    
    (hs[17])->Add(hPlog[0]);
    (hs[17])->Add(hPlog[4]);

    (hs[18])->Add(hPlog[3]);
    (hs[18])->Add(hPlog[5]);

    //PLOT IN TWO CANVAS
    
 
    /*
    for (int i=0, j=0; i<2; i++,j++) {
        (d[i])->Divide(2,2);
        (d[i])->cd(1);
        (hRi[j])->Draw(); //up left ele event
        (d[i])->cd(2);
        (hRi[j+2])->Draw(); //up right muon event
        (d[i])->cd(3);
        (hRi[j+4])->Draw(); //down left tau event
        
        (d[i])->cd(4);
        (hs[j+7])->Draw(); //down right stacked
    }
    */
    
    d[0]->cd();
    hs[7]->Draw();
    hs[7]->GetXaxis()->SetTitle("p_{T} [GeV]");
    hs[7]->GetYaxis()->SetTitle("Events");

    auto legend_ele = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend_ele->AddEntry(hRi[0],"Electron Events","f");
    legend_ele->AddEntry(hRi[2],"Muon Events","f");
    legend_ele->AddEntry(hRi[4],"Tau Events","f");
    legend_ele->Draw();

    d[0]->Modified();
    d[0]->SaveAs("./Graphs/new/Electron_pT_stacked_RctIso.pdf");
    
    
    
    
    
    
    d[1]->cd();
    hs[8]->Draw();
    hs[8]->GetXaxis()->SetTitle("p_{T} [GeV]");
    hs[8]->GetYaxis()->SetTitle("Events/GeV");
    
    auto legend_muon = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend_muon->AddEntry(hRi[1],"Electron Events","f");
    legend_muon->AddEntry(hRi[3],"Muon Events","f");
    legend_muon->AddEntry(hRi[5],"Tau Events","f");
    
    legend_muon->Draw();

    d[1]->Modified();
    d[1]->SaveAs("./Graphs/new/Muon_pT_stacked_RctIso.pdf");
    
    
    
    
    d[2]->cd();
    hs[9]->Draw();
    hs[9]->GetXaxis()->SetTitle("dxy wrt first PV [cm]");
    hs[9]->GetYaxis()->SetTitle("Events");
    
    auto legend_pe = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend_pe->AddEntry(hP[0],"Electron Events","f");
    //legend_ele->AddEntry(hP[2],"Muon Events","f");
    legend_pe->AddEntry(hP[4],"Tau Events","f");
    legend_pe->Draw();
    
    d[2]->SetLogy();
    d[2]->Modified();
    d[2]->SaveAs("./Graphs/new/Electron_dxy.png");
    d[2]->SaveAs("./Graphs/new/Electron_dxy.pdf");
    
    
    
    d[3]->cd();
    hs[10]->Draw();
    hs[10]->GetXaxis()->SetTitle("dxy wrt first PV [cm]");
    hs[10]->GetYaxis()->SetTitle("Events");
    
    auto legend_pm = new TLegend(0.6, 0.7, 0.9, 0.9);
    //legend_ele->AddEntry(hP[0],"Electron Events","f");
    legend_pm->AddEntry(hP[3],"Muon Events","f");
    legend_pm->AddEntry(hP[5],"Tau Events","f");
    legend_pm->Draw();
    
    d[3]->SetLogy();
    d[3]->Modified();
    d[3]->SaveAs("./Graphs/new/Muon_dxy.png");
    d[3]->SaveAs("./Graphs/new/Muon_dxy.pdf");
    
    
    
    d[4]->cd();
    hs[11]->Draw();
    hs[11]->GetXaxis()->SetTitle("dxy/dxyErr wrt first PV [cm]");
    hs[11]->GetYaxis()->SetTitle("Events");
    
    auto legend_peerr = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend_peerr->AddEntry(hPerr[0],"Electron Events","f");
    //legend_ele->AddEntry(hP[2],"Muon Events","f");
    legend_peerr->AddEntry(hPerr[4],"Tau Events","f");
    legend_peerr->Draw();
    
    d[4]->SetLogy();
    d[4]->Modified();
    d[4]->SaveAs("./Graphs/new/Electron_dxyErr.png");
    d[4]->SaveAs("./Graphs/new/Electron_dxyErr.pdf");
    
    
    
    d[5]->cd();
    hs[12]->Draw();
    hs[12]->GetXaxis()->SetTitle("dxy/dxyErr wrt first PV [cm]");
    hs[12]->GetYaxis()->SetTitle("Events");
    
    auto legend_pmerr = new TLegend(0.6, 0.7, 0.9, 0.9);
    //legend_ele->AddEntry(hP[0],"Electron Events","f");
    legend_pmerr->AddEntry(hPerr[3],"Muon Events","f");
    legend_pmerr->AddEntry(hPerr[5],"Tau Events","f");
    legend_pmerr->Draw();

    d[5]->SetLogy();
    d[5]->Modified();
    d[5]->SaveAs("./Graphs/new/Muon_dxyErr.png");
    d[5]->SaveAs("./Graphs/new/Muon_dxyErr.pdf");

    
    
    d[6]->cd();
    hs[13]->Draw();
    hs[13]->GetXaxis()->SetTitle("dz wrt first PV [cm]");
    hs[13]->GetYaxis()->SetTitle("Events");
    
    auto legend_pze = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend_pze->AddEntry(hPz[0],"Electron Events","f");
    //legend_ele->AddEntry(hP[2],"Muon Events","f");
    legend_pze->AddEntry(hPz[4],"Tau Events","f");
    legend_pze->Draw();
    
    d[6]->SetLogy();
    d[6]->Modified();
    d[6]->SaveAs("./Graphs/new/Electron_dz.png");
    d[6]->SaveAs("./Graphs/new/Electron_dz.pdf");
    
    
    
    d[7]->cd();
    hs[14]->Draw();
    hs[14]->GetXaxis()->SetTitle("dz wrt first PV [cm]");
    hs[14]->GetYaxis()->SetTitle("Events");
    
    auto legend_pzm = new TLegend(0.6, 0.7, 0.9, 0.9);
    //legend_ele->AddEntry(hP[0],"Electron Events","f");
    legend_pzm->AddEntry(hPz[3],"Muon Events","f");
    legend_pzm->AddEntry(hPz[5],"Tau Events","f");
    legend_pzm->Draw();
    
    d[7]->SetLogy();
    d[7]->Modified();
    d[7]->SaveAs("./Graphs/new/Muon_dz.png");
    d[7]->SaveAs("./Graphs/new/Muon_dz.pdf");
    
    
    
    d[8]->cd();
    hs[15]->Draw();
    hs[15]->GetXaxis()->SetTitle("dz/dzErr wrt first PV [cm]");
    hs[15]->GetYaxis()->SetTitle("Events");
    
    auto legend_pezerr = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend_pezerr->AddEntry(hPzerr[0],"Electron Events","f");
    //legend_ele->AddEntry(hP[2],"Muon Events","f");
    legend_pezerr->AddEntry(hPzerr[4],"Tau Events","f");
    legend_pezerr->Draw();
    
    d[8]->SetLogy();
    d[8]->Modified();
    d[8]->SaveAs("./Graphs/new/Electron_dzErr.png");
    d[8]->SaveAs("./Graphs/new/Electron_dzErr.pdf");
    
    
    
    d[9]->cd();
    hs[16]->Draw();
    hs[16]->GetXaxis()->SetTitle("dz/dzErr wrt first PV [cm]");
    hs[16]->GetYaxis()->SetTitle("Events");
    
    auto legend_pmzerr = new TLegend(0.6, 0.7, 0.9, 0.9);
    //legend_ele->AddEntry(hP[0],"Electron Events","f");
    legend_pmzerr->AddEntry(hPzerr[3],"Muon Events","f");
    legend_pmzerr->AddEntry(hPzerr[5],"Tau Events","f");
    legend_pmzerr->Draw();
    
    d[9]->SetLogy();
    d[9]->Modified();
    d[9]->SaveAs("./Graphs/new/Muon_dzErr.png");
    d[9]->SaveAs("./Graphs/new/Muon_dzErr.pdf");
    
    
    
    d[10]->cd();
    hs[17]->Draw();
    hs[17]->GetXaxis()->SetTitle("Log(dxy) wrt first PV [cm]");
    hs[17]->GetYaxis()->SetTitle("Events (log binned)");
    
    auto legend_ploge = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend_ploge->AddEntry(hPlog[0],"Electron Events","f");
    //legend_ele->AddEntry(hP[2],"Muon Events","f");
    legend_ploge->AddEntry(hPlog[4],"Tau Events","f");
    legend_ploge->Draw();
    
    d[10]->Modified();
    d[10]->SaveAs("./Graphs/new/Electron_dxylog.png");
    d[10]->SaveAs("./Graphs/new/Electron_dxylog.pdf");
    
    
    
    d[11]->cd();
    hs[18]->Draw();
    hs[18]->GetXaxis()->SetTitle("Log(dxy) wrt first PV");
    hs[18]->GetYaxis()->SetTitle("Events (log binned)");
    
    auto legend_plogm = new TLegend(0.6, 0.7, 0.9, 0.9);
    //legend_ele->AddEntry(hP[0],"Electron Events","f");
    legend_plogm->AddEntry(hPlog[3],"Muon Events","f");
    legend_plogm->AddEntry(hPlog[5],"Tau Events","f");
    legend_plogm->Draw();
    
    d[11]->Modified();
    d[11]->SaveAs("./Graphs/new/Muon_dxylog.png");
    d[11]->SaveAs("./Graphs/new/Muon_dxylog.pdf");
    
    
}

void Events::WriteToFile() //Write to File function
{
    for (int i=0; i<3; i++) { //pt
        ((h[i]))->Write();
    }
    
    for (int i=0; i<6; i++) { //rec prtc
        ((hR[i]))->Write();
    }
    
    for (int i=0; i<6; i++) { //rec iso prtc
        ((hRi[i]))->Write();
    }
    
    for (int i=0; i<9; i++) { //stacked histo
        ((hs[i]))->Write();
    }
    //top->Write(); //write tree
    
    output->Write(); //write file
    output->Close();
    return;
}


void Events::Impact()
{
    int Nmc=214650;
    int sigmattbar=370000;
    int lumidati=138;
    
    double lumimc=0.580135;
    
    //lumimc=(Nmc)/(sigmattbar);
    
    double korr=lumidati/lumimc;
    
    auto hPmm = new TH1F("hPmm","Impact Parameter of Muon, from muon events",10, 0, 0.20);
    auto hPmt = new TH1F("hPmt","Impact Parameter of Muon, from tau events",10, 0, 0.20);
    auto hPem = new TH1F("hPem","Impact Parameter of Muon, from ele events",10, 0, 0.20);

    
    int Nbins=hPmm->GetNbinsX();
    float temp=0;

    for (int i = 1; i<=Nbins; i++) {
        temp = korr * hP[3]->GetBinContent(i);
        hPmm->SetBinContent(i, temp);
        
        temp = korr * hP[5]->GetBinContent(i);
        hPmt->SetBinContent(i, temp);
        
        temp= korr * data->GetBinContent(i);
        data->SetBinContent(i, temp);
        
        temp= korr * hP[1]->GetBinContent(i);
        hPem->SetBinContent(i, temp);
    }
    
    std::cout << " lumi_mc: " <<(double) lumimc<<std::endl;
    std::cout << " lumi_dati: " <<(double) lumidati<<std::endl;
    std::cout << " korr: " <<(double) korr<<std::endl;

    hPmt->SetStats(kFALSE);
    hPmm->SetStats(kFALSE);
    hPem->SetStats(kFALSE);
    


    //RECOLOR THE HISTOGRAM BY LEPT_TYPE
    ((hP[0]))->SetFillColor(38);
    ((hPem))->SetFillColor(38);

    ((hP[2]))->SetFillColor(45);
    ((hPmm))->SetFillColor(45);

    ((hP[4]))->SetFillColor(30);
    ((hPmt))->SetFillColor(30);
    
    //REFILL
    
    
    (hs[9])->Add(hP[0]);
    (hs[9])->Add(hP[4]);

    (hs[10])->Add(hPmm);
    (hs[10])->Add(hPmt);

    
    //int Nbins=10;
    
    double mu_ev=0;
    double tau_ev=0;
    double inv_err=0;
    
    for(int i=0;i<Nbins;i++){
        mu_ev=hPmm->GetBinContent(i);
        tau_ev=hPmt->GetBinContent(i);
        if(((mu_ev+tau_ev)!=0)&&(tau_ev!=0))
            inv_err=inv_err+(std::pow(tau_ev,2)/(mu_ev+tau_ev));

    }
    
    
    std::cout << " Error: " <<(double) std::sqrt(1/inv_err)<<std::endl;
    
    
    
    
    d[0]->cd();
    hPmm->Draw();
    hPmm->GetXaxis()->SetTitle("dxy wrt first PV [cm]");
    hPmm->GetYaxis()->SetTitle("Muon events");
    
    
    d[0]->SetLogy();
    d[0]->Modified();
    d[0]->SaveAs("./Graphs/new/Muon_fMuon_dxy.pdf");
    
    
    
    d[1]->cd();
    hPmt->Draw();
    hPmt->GetXaxis()->SetTitle("dxy wrt first PV [cm]");
    hPmt->GetYaxis()->SetTitle("Tau Events");
    
    d[1]->SetLogy();
    d[1]->Modified();
    d[1]->SaveAs("./Graphs/new/Muon_fTau_dxy.pdf");
    
    

    
    d[3]->cd();
    hs[10]->Draw();
    hs[10]->GetXaxis()->SetTitle("dxy wrt first PV [cm]");
    hs[10]->GetYaxis()->SetTitle("Events");
    
    auto legend_pm = new TLegend(0.6, 0.7, 0.9, 0.9);
    //legend_ele->AddEntry(hP[0],"Electron Events","f");
    legend_pm->AddEntry(hPmm,"Muon Events","f");
    legend_pm->AddEntry(hPmt,"Tau Events","f");
    legend_pm->Draw();
    
    d[3]->SetLogy();
    d[3]->Modified();
    d[3]->SaveAs("./Graphs/new/Muon_new_dxy.pdf");
    
    
    
    TObjArray *mc = new TObjArray(2);        // MC histograms are put in this array
    mc->Add(hPmm);
    mc->Add(hPmt);
    //mc->Add(hPem);
    
    TFractionFitter* fit = new TFractionFitter(data, mc); // initialise
    //fit->Constrain(0,0,1);               // constrain fraction 0 to be between 0 and 1
    //fit->Constrain(1,0,1);               // constrain fraction 1 to be between 0 and 1
    //fit->SetRangeX(1,28);                    // use only the first 15 bins in the fit
    Int_t status = fit->Fit();               // perform the fit
    std::cout << "fit status: " << status << std::endl;
    if (status == 0) {                       // check on fit status
        d[12]->cd();
        TH1F* result = (TH1F*) fit->GetPlot();
        data->Draw("Ep");
        result->Draw("same");
        d[12]->SetLogy();
        d[12]->Modified();
        d[12]->SaveAs("./Graphs/new/Fit_Imapct.pdf");

    }
    
}


void Events::Loop() //main fucntion, actually
{
    
    
    
    if (fChain == 0) return;
    
    Long64_t nentries = fChain->GetEntriesFast();
    

    
    int ilept=-1; // lepton index in LHEPart collection
    


    Long64_t nbytes = 0, nb = 0;
    
    //nentries=100;
    
    for (Long64_t jentry=0; jentry<nentries;jentry++) { //loop over entries
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        int ilept=TypeIndex(); //Get index of particle, if lepton
        if(ilept<0) break; //break loop if not lepton

        int lept_type=0;
        lept_type=abs(LHEPart_pdgId[ilept]); // lepton pdgID type
        
        
        switch (lept_type) {
            case 11: //if electron

                nele++; //increment electron counter
                Filling(0, ilept); //fill histograms and tree
                break;
                
            case 13: //if muon
                
                nmu++;
                Filling(1, ilept);
                break;
                
            case 15: //if tau

                ntau++;
                Filling(2, ilept);
                break;
                
            default:
                
                break;
                
                
        }
        

        nev++; //increment lepton type events counter
        
    }


}


int main() {
    
    Events t;
    
    t.Loop();
    
    //t.ReconStack(); //do the first
    
    
    //t.IsoStack();
    t.Impact();
    
    //t.WriteToFile();
    
    
    std::cout <<" TOTALE EVENTI " <<nev<<std::endl;
    std::cout << " Numero muoni " <<nmu<< " | n_mu/n_tau="<< (double) nmu/ntau<<std::endl;
    std::cout << " Numero elettroni " <<nele<< " | n_ele/n_tau="<< (double) nele/ntau<<std::endl;
    std::cout << " Numero tau " <<ntau<< " | 2*n_tau/(n_ele+n_mu)="<< (double) 2*ntau/(nele+nmu)<<std::endl;
}


