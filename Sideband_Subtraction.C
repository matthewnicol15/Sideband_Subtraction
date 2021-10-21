
{

  // Open file and get the mass distribution histogram
  TFile *f = new TFile("/mnt/g/Strangeness_3_Topo_3_RGB_Inc_Pass1_Inbending_e_3Kp_No_Low_Moment_Sigma_235_ssd_181220_01.root");
  TH1F *hmass = (TH1F*)f->Get("hmass_kp_3"); // Calculated mass of kaon
  TH1F *hmiss_5_sig = (TH1F*)f->Get("hmiss_5_sig"); // Missing mass plot for "signal" region of kaon mass
  TH1F *hmiss_5_back = (TH1F*)f->Get("hmiss_5_back"); // Missing mass plot for "backgroud" region of kaon mass

  // Rebin the mass histogram
  hmass->Rebin(4);

  // Creating histograms for results plots
  auto *hmiss_5_result = new TH1F("hmiss_5_result","Sideband Subtracted",200,0,4);
  hmiss_5_result->GetXaxis()->SetLabelSize(0.05);
  hmiss_5_result->GetXaxis()->SetTitleSize(0.9);
  hmiss_5_result->GetYaxis()->SetLabelSize(0.05);
  hmiss_5_result->GetYaxis()->SetTitleSize(0.05);

  auto *hmiss_5_result_2 = new TH1F("hmiss_5_result_2","Sideband Subtracted",200,0,4);
  hmiss_5_result_2->GetXaxis()->SetLabelSize(0.05);
  hmiss_5_result_2->GetXaxis()->SetTitleSize(0.9);
  hmiss_5_result_2->GetYaxis()->SetLabelSize(0.05);
  hmiss_5_result_2->GetYaxis()->SetTitleSize(0.05);

  // Defining functions to fit kaon mass plot
  TF1 *f1 = new TF1("f1","pol3(0) + gaus(4)",0.365,0.65); // Total function
  TF1 *f2 = new TF1("f2","gaus(0)",0.365,0.65); // Signal function
  TF1 *f3 = new TF1("f3","pol3(0)",0.365,0.65); // background function

  // Setting some parameters to better fit functions
  // f1->SetParameter(0,-55000);
  // f1->SetParameter(1,550000);
  // f1->SetParameter(2,-1000000);
  // f1->SetParameter(3,640000);
  f1->SetParameter(4,50); // Gaussian height
  f1->FixParameter(5,0.497); // Gaussian mean (kaon mass)
  f1->SetParameter(6,0.02); // Gaussian sigma



  hmass->SetMinimum(0);
  hmass->Draw();

  f1->SetLineColor(kRed);

  // Fitting the total function to the kaon mass plot
  hmass->Fit("f1","RB");

  // Fixing the parameters of the signal function to those from the total fit
  f2->FixParameter(0,f1->GetParameter(4));
  f2->FixParameter(1,f1->GetParameter(5));
  f2->FixParameter(2,f1->GetParameter(6));

  // Fixing the parameters of the background function to those from the total fit
  f3->FixParameter(0,f1->GetParameter(0));
  f3->FixParameter(1,f1->GetParameter(1));
  f3->FixParameter(2,f1->GetParameter(2));
  f3->FixParameter(3,f1->GetParameter(3));
  f2->SetLineColor(kBlack);
  f3->SetLineColor(kBlue);

  // Drawing the signal and backgroud functions on top of mass plot
  f2->Draw("same");
  f3->Draw("same");

  // Defning the different limits and parameters need to set the ranges and calculate the integrals
  Double_t Sigma = fabs(f1->GetParameter(6));
  Double_t Mean = f1->GetParameter(5);
  Double_t Below_Lower_Limit, Below_Upper_Limit, Above_Lower_Limit, Above_Upper_Limit,Signal_Lower_Limit, Signal_Upper_Limit;
  Below_Lower_Limit = Mean - (8 * Sigma);
  Below_Upper_Limit = Mean - (4 * Sigma);
  Above_Lower_Limit = Mean + (4 * Sigma);
  Above_Upper_Limit = Mean + (8 * Sigma);
  Signal_Lower_Limit = Mean - (2 * Sigma);
  Signal_Upper_Limit = Mean + (2 * Sigma);
  cout<<Below_Lower_Limit<<"  "<<Below_Upper_Limit<<"  "<<Above_Lower_Limit<<"  "<<Above_Upper_Limit<<"  "<<Signal_Lower_Limit<<"  "<<Signal_Upper_Limit<<endl;

  // Determining the size of signal and background
  Double_t Background_Below = f3->Integral(Below_Lower_Limit,Below_Upper_Limit);
  Double_t Background_Above = f3->Integral(Above_Lower_Limit,Above_Upper_Limit);
  Double_t Signal_back = f3->Integral(Signal_Lower_Limit,Signal_Upper_Limit);
  Double_t Signal = f2->Integral(Signal_Lower_Limit,Signal_Upper_Limit);
  Double_t Sidebands = Background_Below + Background_Above;
  Double_t Overall = f1->Integral(Signal_Lower_Limit,Signal_Upper_Limit);
  Double_t Ratio_1 = Signal_back / Overall;
  Double_t Ratio_2 = hmiss_5_back->Integral() / hmiss_5_sig->Integral();

  cout<<Ratio_1<<" "<<Ratio_2<<endl;
  cout<<Sidebands<<" "<<Signal<<" "<<Signal_back<<endl;

  // hmiss_5_back->Scale(Ratio_1/Ratio_2);
  Double_t Ratio_3 = hmiss_5_back->Integral() / hmiss_5_sig->Integral();
  cout<<Ratio_3<<endl;

  // Rebinning just to improve fit
  *hmiss_5_sig->Rebin(4);
  *hmiss_5_back->Rebin(4);
  *hmiss_5_result->Rebin(4);
  *hmiss_5_result_2->Rebin(4);
  // hmiss_5_back->Smooth(1000);

  // hmiss_5_sig->SetLineColor(kRed);
  // hmiss_5_sig->Draw();
  // hmiss_5_back->Draw("same");
  TF1 *f5=new TF1("f5","pol1(0)",0.5,4);
  hmiss_5_back->Fit("f5","RB");
  *hmiss_5_result = *hmiss_5_sig - *hmiss_5_back;

  //
  for(Int_t j=1;j<hmiss_5_result_2->GetNbinsX()+1;j++){
    if(hmiss_5_sig->GetBinContent(j)>0)hmiss_5_result_2->SetBinContent(j,hmiss_5_sig->GetBinContent(j) - f5->Eval(hmiss_5_sig->GetBinCenter(j)));
  }

  hmiss_5_result->GetXaxis()->SetNdivisions(505);
  hmiss_5_result->GetYaxis()->SetNdivisions(505);
  hmiss_5_result->SetTitle("MM(e' K^{+} K^{+} K^{+});MM(e' K^{+} K^{+} K^{+}) [GeV];Counts");
  hmiss_5_result->GetXaxis()->SetTitleSize(0.06);
  hmiss_5_result->GetYaxis()->SetTitleSize(0.06);
  hmiss_5_result->GetXaxis()->SetLabelSize(0.06);
  hmiss_5_result->GetYaxis()->SetLabelSize(0.06);

  // Creating lines highlighting particle rest masses
  auto* l1=new TLine(0,0,4,0);
  auto* l2=new TLine(1.322,-70,1.322,270);
  auto* l3=new TLine(1.535,-70,1.535,270);
  l1->SetLineColor(kBlue);
  l2->SetLineColor(kGray);
  l3->SetLineColor(kGray);
  l2->SetLineWidth(3);
  l2->SetLineStyle(7);
  l3->SetLineWidth(3);
  l3->SetLineStyle(7);

  auto *c1=new TCanvas("c1","Missing Mass",1200,800);
  gStyle->SetTitleSize(0.03);
  hmiss_5_result_2->SetLineColor(kRed);
  hmiss_5_result_2->SetMarkerColor(kRed);
  // hmiss_5_result->Rebin(10);
  // hmiss_5_result_2->Rebin(10);
  // hmiss_5_result_2->Scale(0.25);
  hmiss_5_result_2->SetMarkerStyle(8);
  hmiss_5_result_2->SetMarkerSize(0.9);
  hmiss_5_result->SetMarkerStyle(8);
  hmiss_5_result->SetMarkerSize(0.9);
  // hmiss_5_sig->Draw("e");
  // hmiss_5_back->Draw("c,hist,same");
  // hmiss_5_result_2->Smooth(5);
  hmiss_5_result->Draw("e");
  hmiss_5_result_2->Draw("c,hist,same");
  l1->Draw("same");
  // l2->Draw("same");
  // l3->Draw("same");
  gStyle->SetOptStat(0);

  auto *c2=new TCanvas("c2","Mass",1200,800);
  hmass->SetMinimum(0);
  hmass->Draw();

  f1->SetLineColor(kRed);

  hmass->Fit("f1","RB");
  f2->Draw("same");
  f3->Draw("same");
  auto *c3=new TCanvas("c3","Mass",1200,800);
  hmiss_5_sig->SetLineColor(kRed);
  hmiss_5_sig->Draw();
  hmiss_5_back->Draw("same");


}
