// ============================================================
//  OutputWriter.cc
//  Pion Cloud Model -- ASCII + (optional) ROOT output
// ============================================================

#include "OutputWriter.hh"
#include "PhysicsParams.hh"
#include "Results.hh"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <ctime>
#include <string>

#ifdef HAVE_ROOT
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TNamed.h"
#endif

namespace PionCloud {

// ============================================================
//  Internal helpers: enum -> human-readable string
// ============================================================
static std::string ffName(FormFactorType t) {
    switch (t) {
    case FF_MONOPOLE:     return "Monopole (typ=0)";
    case FF_DIPOLE:       return "Dipole (typ=1)";
    case FF_EXP_S:        return "Exponential s-dep (typ=2)";
    case FF_COV_DIPOLE:   return "Covariant Dipole (typ=3)";
    case FF_DIPOLE_S_CH:  return "Dipole s-channel (typ=4)";
    case FF_EXP_T:        return "Exponential t-dep (typ=5)";
    case FF_PAULI_VILLAR: return "Pauli-Villar (typ=6)";
    default:              return "Unknown";
    }
}
static std::string flagName(ContributionFlag f) {
    switch (f) {
    case FLAG_PION:      return "pi-N (flag=0)";
    case FLAG_RHO:       return "rho-N (flag=1)";
    case FLAG_PI_DELTA:  return "pi-Delta (flag=2)";
    default:             return "Unknown";
    }
}
static std::string disName(DisChannel d) {
    switch (d) {
    case DIS_CHARGE_EXCHANGE: return "Charge exchange (dis=0)";
    case DIS_NEUTRAL:         return "Neutral (dis=1)";
    default:                  return "Unknown";
    }
}
static std::string outModeName(int m) {
    switch (m) {
    case OUT_NONE:  return "none";
    case OUT_ASCII: return "ascii";
    case OUT_ROOT:  return "root";
    case OUT_ALL:   return "all";
    default:        return std::to_string(m);
    }
}

// ----------------------------------------------------------
//  channelTag  -- plain ASCII name used for ROOT object names
//                 and file stems (no special characters)
//  channelLabel -- ROOT TLatex markup for display in plots
// ----------------------------------------------------------
static std::string channelTag(ContributionFlag flag, DisChannel dis) {
    std::string ch;
    switch (flag) {
    case FLAG_PION:      ch = "piN";      break;
    case FLAG_RHO:       ch = "rhoN";     break;
    case FLAG_PI_DELTA:  ch = "piDelta";  break;
    default:             ch = "unknown";  break;
    }
    if (dis == DIS_NEUTRAL) ch += "_neutral";
    return ch;
}
// channelLabel is only used in writeROOT (requires HAVE_ROOT)
#ifdef HAVE_ROOT
static std::string channelLabel(ContributionFlag flag, DisChannel dis) {
    std::string ch;
    switch (flag) {
    case FLAG_PION:      ch = "#piN";       break;
    case FLAG_RHO:       ch = "#rhoN";      break;
    case FLAG_PI_DELTA:  ch = "#pi#Delta";  break;
    default:             ch = "unknown";    break;
    }
    if (dis == DIS_NEUTRAL) ch += "^{0}";
    return ch;
}
#endif

// ============================================================
//  Constructor
// ============================================================
OutputWriter::OutputWriter(const RunParams &p, const std::string &base)
    : p_(p), base_(base) {}

// ============================================================
//  write  -- dispatch to ASCII and/or ROOT
// ============================================================
void OutputWriter::write(const ScanResults &results) const
{
    if (p_.outputMode == OUT_NONE) {
        std::cout << "[OutputWriter] Output mode = none. No files written.\n";
        return;
    }
    if (p_.outputMode & OUT_ASCII)
        writeASCII(results);
    if (p_.outputMode & OUT_ROOT)
        writeROOT(results);
}

// ============================================================
//  writeASCII
//
//  Files written:
//    <base>_params.txt     full parameter record + first moments
//    <base>_xscan.dat      x-scan data table
//    <base>_thetascan.dat  theta-scan data table
// ============================================================
void OutputWriter::writeASCII(const ScanResults &results) const
{
    // ---- parameter / summary file ----
    {
        std::string fname = base_ + "_params.txt";
        std::ofstream f(fname);
        if (!f) { std::cerr << "[ASCII] Cannot open " << fname << "\n"; return; }

        std::time_t now = std::time(nullptr);
        char tbuf[64];
        std::strftime(tbuf, sizeof(tbuf), "%Y-%m-%d %H:%M:%S", std::localtime(&now));

        f << "# ============================================================\n"
          << "# Pion Cloud Model -- Run Parameter Record\n"
          << "# Generated : " << tbuf << "\n"
          << "# ============================================================\n\n"

          << "[Run Control]\n"
          << "  scan_mode       = " << p_.scanMode    << "\n"
          << "  output_mode     = " << outModeName(p_.outputMode) << "\n"
          << "  output_name     = " << p_.outBaseName << "\n\n"

          << "[Beam Kinematics]\n"
          << "  beam_energy     = " << p_.E        << " GeV\n"
          << "  electron_angle  = " << p_.theta_e  << " deg\n\n"

          << "[Hadron Kinematics]\n"
          << "  hadron_momentum  = " << p_.pH         << " GeV\n"
          << "  hadron_angle_min = " << p_.alpha1_deg << " deg\n"
          << "  hadron_angle_max = " << p_.alpha2_deg << " deg\n\n"

          << "[Physics Parameters]\n"
          << "  contribution    = " << flagName(p_.flag) << "\n"
          << "  dissociation    = " << disName(p_.dis)   << "\n"
          << "  form_factor     = " << ffName(p_.typ)    << "\n"
          << "  lambda          = " << p_.L   << " GeV  (pi/rho-N vertex)\n"
          << "  lambda_delta    = " << p_.Ld  << " GeV  (Delta vertex)\n\n"

          << "[x-Scan Settings]\n"
          << "  x_min              = " << p_.xmin    << "\n"
          << "  x_max              = " << p_.xmax    << "\n"
          << "  nx                 = " << p_.nx      << "  (output points)\n"
          << "  n_theta_integration= " << p_.nth_int << "  (internal theta steps)\n"
          << "  theta_min          = " << p_.thmin_deg << " deg\n"
          << "  theta_max          = " << p_.thmax_deg << " deg\n"
          << "  k_min              = " << p_.kmin << " GeV  (BoNuS tagging window)\n"
          << "  k_max              = " << p_.kmax << " GeV\n"
          << "  nk                 = " << p_.nk   << "  (integration steps)\n\n"

          << "[theta-Scan Settings]\n"
          << "  theta_min          = " << p_.thmin_deg << " deg\n"
          << "  theta_max          = " << p_.thmax_deg << " deg\n"
          << "  nth                = " << p_.nth    << "  (output points)\n"
          << "  n_x_integration    = " << p_.nx_int << "  (internal x steps)\n"
          << "  x range (internal) = [" << p_.xmin_th << ", " << p_.xmax_th << "]\n"
          << "  k range (internal) = [" << p_.kmin_th << ", " << p_.kmax_th << "] GeV\n"
          << "  nk                 = " << p_.nk << "  (shared with x-scan)\n\n"

          << "[Results]\n"
          << std::scientific << std::setprecision(8)
          << "  First moment F2piK (x-scan)     = " << results.F2pi0_x     << "\n"
          << "  First moment F2piK (theta-scan) = " << results.F2pi0_theta << "\n"
          << "\n[CTEQ6 Neutron F2]\n"
          << "  Available = "
#ifdef HAVE_CTEQ6
          << "yes (CTEQ6L, charge symmetry, no strange)\n";
#else
          << "no (rebuild with -DUSE_CTEQ6=ON to enable F2n and ratio)\n";
#endif

        std::cout << "[ASCII] Written " << fname << "\n";
    }

    // ---- x-scan data file ----
    if (!results.xScan.empty()) {
        // File name encodes channel and mode
        std::string pfxAsc = flagName(p_.flag) + "_" + p_.scanMode;
        std::string fname = base_ + "_" + pfxAsc + ".dat";
        std::ofstream f(fname);
        if (!f) { std::cerr << "[ASCII] Cannot open " << fname << "\n"; return; }

        bool hasCTEQ6 = false;
        for (auto &pt : results.xScan)
            if (pt.F2n > 0.0) { hasCTEQ6 = true; break; }

        std::string nucleonLabel = (p_.nucleon == NUCLEON_PROTON) ? "proton" : "neutron";

        f << "# PionCloud -- x-scan\n"
          << "# channel=" << flagName(p_.flag)
          << "  mode=" << p_.scanMode
          << "  FF=" << ffName(p_.typ)
          << "  L=" << p_.L << " GeV"
          << "  dis=" << disName(p_.dis) << "\n";
        if (p_.scanMode == "kbin")
            f << "# k=[" << p_.kmin << "," << p_.kmax << "] GeV"
              << "  cosph=[" << p_.cosph_min << "," << p_.cosph_max << "]\n";
        if (hasCTEQ6)
            f << "# F2N: CTEQ6M " << nucleonLabel << " F2 via charge symmetry\n"
              << "# ratio: F2pi_" << flagName(p_.flag) << " / F2_" << nucleonLabel << "\n";
        else
            f << "# F2n/ratio: not available (built without CTEQ6)\n";
        f << "#\n"
          << std::setw(14) << "x"
          << std::setw(16) << "Q2_GeV2"
          << std::setw(22) << ("F2piK_" + pfxAsc)
          << std::setw(20) << "RunningMoment";
        if (hasCTEQ6)
            f << std::setw(14) << ("F2_" + nucleonLabel)
              << std::setw(20) << ("ratio_" + pfxAsc);
        f << "\n"
          << std::scientific << std::setprecision(8);
        for (auto &pt : results.xScan) {
            f << std::setw(14) << pt.x
              << std::setw(16) << pt.Q2
              << std::setw(22) << pt.F2piK
              << std::setw(20) << pt.F2piMoment;
            if (hasCTEQ6)
                f << std::setw(14) << pt.F2n
                  << std::setw(20) << pt.ratio;
            f << "\n";
        }
        std::cout << "[ASCII] Written " << fname << "\n";
    }

    // ---- theta-scan data file ----
    if (!results.thetaScan.empty()) {
        std::string fname = base_ + "_thetascan.dat";
        std::ofstream f(fname);
        if (!f) { std::cerr << "[ASCII] Cannot open " << fname << "\n"; return; }

        std::string tpfxAsc = channelTag(p_.flag, p_.dis) + "_theta";
        f << "# PionCloud -- theta_h-scan\n"
          << "# channel=" << flagName(p_.flag)
          << "  dis=" << disName(p_.dis)
          << "  FF=" << ffName(p_.typ)
          << "  L=" << p_.L << " GeV\n"
          << "# x in [" << p_.xmin_th << "," << p_.xmax_th << "]"
          << "  k in [" << p_.kmin_th << "," << p_.kmax_th << "] GeV\n"
          << "#\n"
          << std::setw(16) << "theta_h_deg"
          << std::setw(26) << ("F2piK_" + tpfxAsc)
          << std::setw(20) << "RunningMoment\n";
        f << std::scientific << std::setprecision(8);
        for (auto &pt : results.thetaScan) {
            f << std::setw(16) << pt.theta_deg
              << std::setw(26) << pt.F2piK
              << std::setw(20) << pt.F2piMoment << "\n";
        }
        std::cout << "[ASCII] Written " << fname << "\n";
    }
}

// ============================================================
//  writeROOT
//
//  All object names encode the channel and scan mode so that
//  multiple ROOT files can be merged without name collisions.
//  The prefix is: <channel>_<mode>
//    channel = piN | piN_neutral | rhoN | piDelta
//    mode    = xscan | theta | kbin
//
//  Example: piN_kbin, rhoN_kbin, piDelta_kbin, piN_neutral_theta
//
//  Objects always written (xScan non-empty):
//    g_<pfx>         TGraph   F₂^channel(x)
//    g_moment_<pfx>  TGraph   running ∫F₂ dx
//    h_<pfx>         TH1D     same as histogram
//    c_<pfx>         TCanvas  single-curve plot (log scale)
//
//  Additional objects when CTEQ6 data is present:
//    g_F2neutron           TGraph  F₂ⁿ(x) -- always same name
//    g_ratio_<pfx>         TGraph  F₂^channel / F₂ⁿ
//    h_ratio_<pfx>         TH1D    ratio as histogram
//    c_ratio_<pfx>         TCanvas two-pad: curves top, ratio bottom
//
//  Objects written when thetaScan non-empty:
//    g_thetascan     TGraph   F₂^π(θ_h)
//    g_moment_theta  TGraph   running ∫F₂ dθ
//    h_thetascan     TH1D     same as histogram
//    c_thetascan     TCanvas  single-curve plot
//
//  RunParams  TNamed  full parameter string
// ============================================================
void OutputWriter::writeROOT(const ScanResults &results) const
{
#ifndef HAVE_ROOT
    (void)results;
    std::cerr << "[ROOT] ROOT output requested but HAVE_ROOT not defined at "
                 "compile time.\n"
              << "       Rebuild with: cmake -DUSE_ROOT=ON  (requires ROOT >= 6)\n"
              << "       or:  g++ ... $(root-config --cflags --libs) -DHAVE_ROOT\n";
#else
    std::string fname = base_ + ".root";
    TFile *rootFile = TFile::Open(fname.c_str(), "RECREATE");
    if (!rootFile || rootFile->IsZombie()) {
        std::cerr << "[ROOT] Cannot create " << fname << "\n";
        return;
    }
    rootFile->cd();

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // ---- name prefix: <channel>_<mode>, e.g. "piN_kbin" ----
    // channelTag()  -> plain ASCII for ROOT object names
    // channelLabel() -> ROOT TLatex for axis titles and legends
    const std::string pfx   = channelTag  (p_.flag, p_.dis) + "_" + p_.scanMode;
    const std::string cTag  = channelTag  (p_.flag, p_.dis);
    const std::string cLbl  = channelLabel(p_.flag, p_.dis);

    // ---- x-scan / kbin-scan objects (both use xScan vector) ----
    if (!results.xScan.empty()) {
        int n = static_cast<int>(results.xScan.size());

        std::vector<double> vx(n), vF(n), vM(n), vFn(n), vR(n);
        bool hasCTEQ6 = false;
        for (int i = 0; i < n; ++i) {
            vx[i] = results.xScan[i].x;
            vF[i] = results.xScan[i].F2piK;
            vM[i] = results.xScan[i].F2piMoment;
            vFn[i]= results.xScan[i].F2n;
            vR[i] = results.xScan[i].ratio;
            if (vFn[i] > 0.0) hasCTEQ6 = true;
        }

        // Build Y-axis label using TLatex channel name
        std::string yLabel;
        if (p_.scanMode == "kbin")
            yLabel = "F_{2}^{(" + cLbl + ")}(x; #Deltak, #Delta#theta_{h})";
        else
            yLabel = "F_{2}^{(" + cLbl + ")}(x)";

        // TGraph: F2pi vs x
        std::string gFname = "g_" + pfx;
        auto *gF = new TGraph(n, vx.data(), vF.data());
        gF->SetName(gFname.c_str());
        gF->GetXaxis()->SetTitle("x_{Bj}");
        gF->GetYaxis()->SetTitle(yLabel.c_str());
        gF->SetLineColor(kBlue+1);   gF->SetLineWidth(2);
        gF->SetMarkerColor(kBlue+1); gF->SetMarkerStyle(20); gF->SetMarkerSize(0.8);
        gF->Write();

        // TGraph: running moment vs x
        std::string gMname = "g_moment_" + pfx;
        auto *gM = new TGraph(n, vx.data(), vM.data());
        gM->SetName(gMname.c_str());
        gM->GetXaxis()->SetTitle("x_{Bj}");
        gM->GetYaxis()->SetTitle(("#int F_{2}^{(" + cLbl + ")} dx").c_str());
        gM->SetLineColor(kRed+1); gM->SetLineWidth(2);
        gM->SetMarkerColor(kRed+1); gM->SetMarkerStyle(24); gM->SetMarkerSize(0.8);
        gM->Write();

        // TH1D
        double halfStep = (n > 1) ? 0.5*(vx[1]-vx[0]) : 0.01;
        std::string hFname = "h_" + pfx;
        auto *hF = new TH1D(hFname.c_str(),
                            (";x_{Bj};" + yLabel).c_str(),
                            n, vx.front()-halfStep, vx.back()+halfStep);
        for (int i = 0; i < n; ++i) hF->SetBinContent(i+1, vF[i]);
        hF->SetLineColor(kBlue+1); hF->SetLineWidth(2);
        hF->Write();

        // CTEQ6-dependent objects
        if (hasCTEQ6) {
            auto *gFn = new TGraph(n, vx.data(), vFn.data());
            gFn->SetName("g_F2neutron");
            gFn->GetXaxis()->SetTitle("x_{Bj}");
            gFn->GetYaxis()->SetTitle("F_{2}^{n}(x)  [CTEQ6M]");
            gFn->SetLineColor(kOrange+2); gFn->SetLineWidth(2);
            gFn->SetMarkerColor(kOrange+2); gFn->SetMarkerStyle(22); gFn->SetMarkerSize(0.8);
            gFn->Write();

            std::string gRname = "g_ratio_" + pfx;
            auto *gR = new TGraph(n, vx.data(), vR.data());
            gR->SetName(gRname.c_str());
            gR->GetXaxis()->SetTitle("x_{Bj}");
            gR->GetYaxis()->SetTitle(("F_{2}^{(" + cLbl + ")} / F_{2}^{n}").c_str());
            gR->SetLineColor(kMagenta+1); gR->SetLineWidth(2);
            gR->SetMarkerColor(kMagenta+1); gR->SetMarkerStyle(21); gR->SetMarkerSize(0.8);
            gR->Write();

            std::string hRname = "h_ratio_" + pfx;
            auto *hR = new TH1D(hRname.c_str(),
                                (";x_{Bj};F_{2}^{(" + cLbl + ")}/F_{2}^{n}").c_str(),
                                n, vx.front()-halfStep, vx.back()+halfStep);
            for (int i = 0; i < n; ++i) hR->SetBinContent(i+1, vR[i]);
            hR->SetLineColor(kMagenta+1); hR->SetLineWidth(2);
            hR->Write();

            std::string cRname = "c_ratio_" + pfx;
            auto *cRatio = new TCanvas(cRname.c_str(),
                                       ("F2pi/F2n: " + pfx).c_str(), 900, 800);
            cRatio->Divide(1, 2, 0.0, 0.0);
            cRatio->cd(1);
            gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.05);
            gPad->SetBottomMargin(0.02); gPad->SetTopMargin(0.10);
            gPad->SetPad(0.0, 0.35, 1.0, 1.0);
            gPad->SetLogy(1);
            gF->Draw("ALP");
            gF->GetXaxis()->SetLabelSize(0);
            gF->GetXaxis()->SetTitleSize(0);
            gF->GetYaxis()->SetTitleSize(0.06);
            gF->GetYaxis()->SetLabelSize(0.05);
            gF->GetYaxis()->SetTitleOffset(1.1);
            gFn->Draw("LP SAME");
            auto *leg = new TLegend(0.55, 0.65, 0.92, 0.88);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.05);
            leg->AddEntry(gF,  yLabel.c_str(), "lp");
            leg->AddEntry(gFn, "F_{2}^{n} (CTEQ6M)", "lp");
            leg->Draw();
            auto *tex1 = new TLatex();
            tex1->SetNDC(); tex1->SetTextSize(0.040); tex1->SetTextColor(kGray+2);
            std::ostringstream lbl1;
            lbl1 << cLbl << ",  " << ffName(p_.typ)
                 << ",  #Lambda=" << p_.L << " GeV";
            if (p_.scanMode == "kbin")
                lbl1 << ",  |k|=[" << p_.kmin*1000 << "," << p_.kmax*1000 << "] MeV";
            tex1->DrawLatex(0.15, 0.92, lbl1.str().c_str());
            cRatio->cd(2);
            gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.05);
            gPad->SetBottomMargin(0.18); gPad->SetTopMargin(0.02);
            gPad->SetPad(0.0, 0.0, 1.0, 0.35);
            gR->Draw("ALP");
            gR->GetXaxis()->SetTitleSize(0.10);
            gR->GetXaxis()->SetLabelSize(0.08);
            gR->GetYaxis()->SetTitleSize(0.09);
            gR->GetYaxis()->SetLabelSize(0.08);
            gR->GetYaxis()->SetTitleOffset(0.65);
            gR->GetYaxis()->SetNdivisions(505);
            cRatio->Write();
            delete cRatio; delete leg; delete tex1;
        }

        std::string cFname = "c_" + pfx;
        auto *c = new TCanvas(cFname.c_str(), pfx.c_str(), 900, 650);
        c->SetLeftMargin(0.14); c->SetRightMargin(0.05);
        c->SetBottomMargin(0.13); c->SetTopMargin(0.08);
        c->SetLogy(1);
        gF->Draw("ALP");
        gF->GetXaxis()->SetTitleSize(0.05);
        gF->GetYaxis()->SetTitleSize(0.05);
        gF->GetXaxis()->SetLabelSize(0.04);
        gF->GetYaxis()->SetLabelSize(0.04);
        gF->GetYaxis()->SetTitleOffset(1.3);
        auto *tex = new TLatex();
        tex->SetNDC(); tex->SetTextSize(0.030); tex->SetTextColor(kGray+2);
        std::ostringstream lbl;
        lbl << cLbl << ",  " << ffName(p_.typ)
            << ",  #Lambda=" << p_.L << " GeV";
        if (p_.scanMode == "kbin")
            lbl << ",  |k|=[" << p_.kmin*1000 << "," << p_.kmax*1000 << "] MeV"
                << ",  #theta=[30,70]#circ";
        tex->DrawLatex(0.15, 0.93, lbl.str().c_str());
        c->Write();
        delete c; delete tex;
    }

    // ---- theta-scan objects ----
    if (!results.thetaScan.empty()) {
        int n = static_cast<int>(results.thetaScan.size());

        // theta-scan prefix follows same convention as x/kbin
        std::string tpfx = cTag + "_theta";
        std::string tyLabel = "F_{2}^{(" + cLbl + ")}(#theta_{h})";

        std::vector<double> vth(n), vF(n), vM(n);
        for (int i = 0; i < n; ++i) {
            vth[i] = results.thetaScan[i].theta_deg;
            vF[i]  = results.thetaScan[i].F2piK;
            vM[i]  = results.thetaScan[i].F2piMoment;
        }

        auto *gF = new TGraph(n, vth.data(), vF.data());
        gF->SetName(("g_" + tpfx).c_str());
        gF->GetXaxis()->SetTitle("#theta_{h} [deg]");
        gF->GetYaxis()->SetTitle(tyLabel.c_str());
        gF->SetLineColor(kGreen+2); gF->SetLineWidth(2);
        gF->SetMarkerColor(kGreen+2); gF->SetMarkerStyle(21); gF->SetMarkerSize(0.8);
        gF->Write();

        auto *gM = new TGraph(n, vth.data(), vM.data());
        gM->SetName(("g_moment_" + tpfx).c_str());
        gM->GetXaxis()->SetTitle("#theta_{h} [deg]");
        gM->GetYaxis()->SetTitle(("#int F_{2}^{(" + cLbl + ")} d#theta").c_str());
        gM->SetLineColor(kOrange+2); gM->SetLineWidth(2);
        gM->SetMarkerColor(kOrange+2); gM->SetMarkerStyle(25); gM->SetMarkerSize(0.8);
        gM->Write();

        double halfStep = (n > 1) ? 0.5*(vth[1]-vth[0]) : 0.5;
        auto *hF = new TH1D(("h_" + tpfx).c_str(),
                            (";#theta_{h} [deg];" + tyLabel).c_str(),
                            n, vth.front()-halfStep, vth.back()+halfStep);
        for (int i = 0; i < n; ++i) hF->SetBinContent(i+1, vF[i]);
        hF->SetLineColor(kGreen+2); hF->SetLineWidth(2);
        hF->Write();

        auto *c = new TCanvas(("c_" + tpfx).c_str(),
                              ("theta-scan: " + tpfx).c_str(), 900, 650);
        c->SetLeftMargin(0.14); c->SetRightMargin(0.05);
        c->SetBottomMargin(0.13); c->SetTopMargin(0.08);
        c->SetLogy(1);
        gF->Draw("ALP");
        gF->GetXaxis()->SetTitleSize(0.05);
        gF->GetYaxis()->SetTitleSize(0.05);
        gF->GetXaxis()->SetLabelSize(0.04);
        gF->GetYaxis()->SetLabelSize(0.04);
        gF->GetYaxis()->SetTitleOffset(1.3);
        auto *tex = new TLatex();
        tex->SetNDC(); tex->SetTextSize(0.030); tex->SetTextColor(kGray+2);
        std::ostringstream lbl;
        lbl << cLbl << ",  " << ffName(p_.typ)
            << ",  #Lambda=" << p_.L << " GeV"
            << ",  #Deltax=[" << p_.xmin_th << "," << p_.xmax_th << "]"
            << ",  #Deltak=[" << p_.kmin_th*1000 << "," << p_.kmax_th*1000 << "] MeV";
        tex->DrawLatex(0.15, 0.93, lbl.str().c_str());
        c->Write();
        delete c; delete tex;
    }

    // ---- Full parameter string stored as TNamed ----
    std::ostringstream ps;
    ps << "scan_mode="     << p_.scanMode
       << "; E="           << p_.E        << " GeV"
       << "; theta_e="     << p_.theta_e  << " deg"
       << "; contribution="<< flagName(p_.flag)
       << "; dis="         << disName(p_.dis)
       << "; FF="          << ffName(p_.typ)
       << "; L="           << p_.L        << " GeV"
       << "; Ld="          << p_.Ld       << " GeV"
       << "; xmin="        << p_.xmin
       << "; xmax="        << p_.xmax
       << "; nx="          << p_.nx
       << "; nth_int="     << p_.nth_int
       << "; kmin="        << p_.kmin     << " GeV"
       << "; kmax="        << p_.kmax     << " GeV"
       << "; nk="          << p_.nk
       << "; thmin="       << p_.thmin_deg<< " deg"
       << "; thmax="       << p_.thmax_deg<< " deg"
       << "; nth="         << p_.nth
       << "; nx_int="      << p_.nx_int
       << "; F2pi0_x="     << results.F2pi0_x
       << "; F2pi0_theta=" << results.F2pi0_theta
#ifdef HAVE_CTEQ6
       << "; CTEQ6=yes"
#else
       << "; CTEQ6=no"
#endif
       ;

    auto *pobj = new TNamed("RunParams", ps.str().c_str());
    pobj->Write();
    delete pobj;

    rootFile->Close();
    delete rootFile;
    std::cout << "[ROOT] Written " << fname << "\n";
#endif
}

} // namespace PionCloud
