/**
 * @file run_all_plots.C
 * @brief A ROOT macro to compare particle acceptance from different data sources.
 * @author Gemini
 * @version 1.7
 * @date 2025-09-17
 *
 * This macro reads TH1F histograms and generates comparison plots.
 * It features an advanced algorithm to find an empty rectangle on the canvas
 * for legend placement, constrained inside the plot axes with a safety buffer.
 * The canvas margins are adjusted for optimal label visibility in PDFs.
 */

#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TPad.h> // Needed to get pad margins
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <map>

/**
 * @brief Finds the quietest quadrant for a legend (Fallback Method).
 */
std::vector<double> findQuietestQuadrant(TH1* h1, TH1* h2) {
    const double x_center = 0.5, y_center = 0.5;
    std::map<std::string, double> scores = {{"TL", 0.0}, {"TR", 0.0}, {"BL", 0.0}, {"BR", 0.0}};
    
    TAxis* xaxis = h1->GetXaxis(); TAxis* yaxis = h1->GetYaxis();
    double x_range = xaxis->GetXmax() - xaxis->GetXmin();
    double y_range = yaxis->GetXmax() - yaxis->GetXmin();

    if (x_range <= 0 || y_range <= 0) return {0.62, 0.65, 0.92, 0.90};

    for (int i = 1; i <= h1->GetNbinsX(); ++i) {
        double x_norm = (xaxis->GetBinCenter(i) - xaxis->GetXmin()) / x_range;
        double y1 = h1->GetBinContent(i) + h1->GetBinError(i);
        double y2 = h2->GetBinContent(i) + h2->GetBinError(i);
        double y_norm = (std::max(y1, y2) - yaxis->GetXmin()) / y_range;
        
        if (x_norm < x_center) {
            if (y_norm > y_center) scores["TL"] += y_norm; else scores["BL"] += y_norm;
        } else {
            if (y_norm > y_center) scores["TR"] += y_norm; else scores["BR"] += y_norm;
        }
    }
    
    double min_score = 1e9;
    std::string best_quadrant = "TR";
    for (auto const& [quadrant, score] : scores) {
        if (score < min_score) {
            min_score = score;
            best_quadrant = quadrant;
        }
    }
    
    // Adjusted coordinates to avoid touching the boundary
    if (best_quadrant == "TL") return {0.20, 0.65, 0.50, 0.90};
    if (best_quadrant == "TR") return {0.62, 0.65, 0.92, 0.90};
    if (best_quadrant == "BL") return {0.20, 0.15, 0.50, 0.40};
    if (best_quadrant == "BR") return {0.62, 0.15, 0.92, 0.40};
    return {0.62, 0.65, 0.92, 0.90};
}

/**
 * @brief Finds an empty, axis-bound rectangle for legend placement.
 */
std::vector<double> findEmptyRectangle(TH1* h1, TH1* h2) {
    const int GRID_X = 25, GRID_Y = 25;
    const double LEGEND_W_NORM = 0.30, LEGEND_H_NORM = 0.20;
    const int LEGEND_W_CELLS = static_cast<int>(LEGEND_W_NORM * GRID_X);
    const int LEGEND_H_CELLS = static_cast<int>(LEGEND_H_NORM * GRID_Y);

    bool occupied[GRID_Y][GRID_X] = {{false}};

    if (!gPad) return findQuietestQuadrant(h1, h2);
    double lm = gPad->GetLeftMargin(), rm = gPad->GetRightMargin();
    double bm = gPad->GetBottomMargin(), tm = gPad->GetTopMargin();
    double safety_buffer = 0.02; // Prevents legend from touching axis frame

    for (int r = 0; r < GRID_Y; ++r) {
        for (int c = 0; c < GRID_X; ++c) {
            double x_norm = (c + 0.5) / GRID_X;
            double y_norm = (r + 0.5) / GRID_Y;
            if (x_norm < lm || x_norm > (1 - rm - safety_buffer) || 
                y_norm < bm || y_norm > (1 - tm - safety_buffer)) {
                occupied[r][c] = true;
            }
        }
    }

    TAxis* xaxis = h1->GetXaxis(); TAxis* yaxis = h1->GetYaxis();
    double x_min_axis = xaxis->GetXmin(), x_max_axis = xaxis->GetXmax();
    double y_min_axis = yaxis->GetXmin(), y_max_axis = yaxis->GetXmax();
    double x_range = x_max_axis - x_min_axis;
    double y_range = y_max_axis - y_min_axis;

    if (x_range <= 0 || y_range <= 0) return findQuietestQuadrant(h1, h2);

    for (TH1* h : {h1, h2}) {
        for (int i = 1; i <= h->GetNbinsX(); ++i) {
            double x_start = xaxis->GetBinLowEdge(i), x_end = xaxis->GetBinUpEdge(i);
            double y_low = h->GetBinContent(i) - h->GetBinError(i);
            double y_high = h->GetBinContent(i) + h->GetBinError(i);

            double x_norm_start = lm + ((x_start - x_min_axis) / x_range) * (1 - lm - rm);
            double x_norm_end = lm + ((x_end - x_min_axis) / x_range) * (1 - lm - rm);
            double y_norm_low = bm + ((y_low - y_min_axis) / y_range) * (1 - bm - tm);
            double y_norm_high = bm + ((y_high - y_min_axis) / y_range) * (1 - bm - tm);

            int gx_start = std::max(0, static_cast<int>(x_norm_start * GRID_X));
            int gx_end = std::min(GRID_X - 1, static_cast<int>(x_norm_end * GRID_X));
            int gy_start = std::max(0, static_cast<int>(y_norm_low * GRID_Y));
            int gy_end = std::min(GRID_Y - 1, static_cast<int>(y_norm_high * GRID_Y));

            for (int r = gy_start; r <= gy_end; ++r) {
                for (int c = gx_start; c <= gx_end; ++c) {
                    occupied[r][c] = true;
                }
            }
        }
    }

    int search_orders[4][2] = {{1,-1}, {-1,-1}, {1,1}, {-1,1}};
    for(int i=0; i<4; ++i) {
        int y_dir = search_orders[i][0], x_dir = search_orders[i][1];
        int r_start = (y_dir > 0) ? 0 : GRID_Y - LEGEND_H_CELLS;
        int r_end = (y_dir > 0) ? GRID_Y - LEGEND_H_CELLS : -1;
        int c_start = (x_dir > 0) ? 0 : GRID_X - LEGEND_W_CELLS;
        int c_end = (x_dir > 0) ? GRID_X - LEGEND_W_CELLS : -1;

        for (int r = r_start; r != r_end; r += y_dir) {
            for (int c = c_start; c != c_end; c += x_dir) {
                bool is_empty = true;
                for (int dr = 0; dr < LEGEND_H_CELLS; ++dr) {
                    for (int dc = 0; dc < LEGEND_W_CELLS; ++dc) {
                        if (occupied[r + dr][c + dc]) {
                            is_empty = false; break;
                        }
                    }
                    if (!is_empty) break;
                }
                if (is_empty) {
                    double x1 = static_cast<double>(c) / GRID_X;
                    double y1 = static_cast<double>(r) / GRID_Y;
                    return {x1, y1, x1 + LEGEND_W_NORM, y1 + LEGEND_H_NORM};
                }
            }
        }
    }
    
    std::cout << "Warning: Could not find a perfectly empty space for the legend in plot '" << h1->GetTitle() << "'. Falling back." << std::endl;
    return findQuietestQuadrant(h1, h2);
}

void plot_LH2_acceptance() {
    TFile *f_mass = TFile::Open("acceptance_mass_xF.root");
    TFile *f_h_shivangi = TFile::Open("acceptance_h_shivangi.root");
    if (!f_mass || f_mass->IsZombie() || !f_h_shivangi || f_h_shivangi->IsZombie()) return;

    for (int i = 0; i < 16; ++i) {
        TCanvas *c = new TCanvas(Form("c_lh2_%d", i), Form("LH2 Acceptance xF bin %d", i), 800, 600);
        c->SetLeftMargin(0.15); // Increase left margin for Y-axis title
        c->SetTickx(1); c->SetTicky(1);
        TH1F *h_mass = (TH1F*)f_mass->Get(Form("h_ratio_LH2_xF_bin%d", i));
        TH1F *h_shivangi = (TH1F*)f_h_shivangi->Get(Form("LH2_acc_M_xF_%d", i));
        if (!h_mass || !h_shivangi) continue;

        h_mass->SetStats(0); h_shivangi->SetStats(0);
        h_mass->SetLineColor(kRed); h_shivangi->SetLineColor(kBlue);
        double maxY = 0;
        for (int bin = 1; bin <= h_mass->GetNbinsX(); ++bin) {
            maxY = std::max({maxY, h_mass->GetBinContent(bin) + h_mass->GetBinError(bin), h_shivangi->GetBinContent(bin) + h_shivangi->GetBinError(bin)});
        }
        h_mass->GetYaxis()->SetRangeUser(0, maxY * 1.1);
        h_mass->SetTitle(Form("LH2 Acceptance for 0.%02d #leq xF < 0.%02d", 5*i, 5*(i+1)));
        h_mass->GetXaxis()->SetTitle("Mass (GeV/c^{2})"); h_mass->GetYaxis()->SetTitle("Acceptance");
        h_mass->GetXaxis()->CenterTitle(true); h_mass->GetYaxis()->CenterTitle(true);
        h_mass->Draw("E1"); h_shivangi->Draw("E1 SAME");
        gPad->Update();

        std::vector<double> leg_coords = findEmptyRectangle(h_mass, h_shivangi);
        //TLegend *leg = new TLegend(leg_coords[0]-0.1, leg_coords[1], leg_coords[2]-0.1, leg_coords[3]);
        TLegend *leg = new TLegend(0.4, leg_coords[1], 0.8, leg_coords[3]);
        leg->SetBorderSize(0);
        leg->AddEntry(h_mass, "Latest Acceptance", "l");
        leg->AddEntry(h_shivangi, "Shivangi's Acceptance", "l");
        leg->Draw();
        c->SaveAs(Form("LH2_acceptance_xF_bin_%d.pdf", i));
        delete c;
    }
    f_mass->Close(); f_h_shivangi->Close();
}

void plot_LD2_acceptance() {
    TFile *f_mass = TFile::Open("acceptance_mass_xF.root");
    TFile *f_d_shivangi = TFile::Open("acceptance_d_shivangi.root");
    if (!f_mass || f_mass->IsZombie() || !f_d_shivangi || f_d_shivangi->IsZombie()) return;

    for (int i = 0; i < 16; ++i) {
        TCanvas *c = new TCanvas(Form("c_ld2_%d", i), Form("LD2 Acceptance xF bin %d", i), 800, 600);
        c->SetLeftMargin(0.15);
        c->SetTickx(1); c->SetTicky(1);
        TH1F *h_mass = (TH1F*)f_mass->Get(Form("h_ratio_LD2_xF_bin%d", i));
        TH1F *h_shivangi = (TH1F*)f_d_shivangi->Get(Form("LD2_acc_M_xF_%d", i));
        if (!h_mass || !h_shivangi) continue;
        
        h_mass->SetStats(0); h_shivangi->SetStats(0);
        h_mass->SetLineColor(kRed); h_shivangi->SetLineColor(kBlue);
        double maxY = 0;
        for (int bin = 1; bin <= h_mass->GetNbinsX(); ++bin) {
            maxY = std::max({maxY, h_mass->GetBinContent(bin) + h_mass->GetBinError(bin), h_shivangi->GetBinContent(bin) + h_shivangi->GetBinError(bin)});
        }
        h_mass->GetYaxis()->SetRangeUser(0, maxY * 1.1);
        h_mass->SetTitle(Form("LD2 Acceptance for 0.%02d #leq xF < 0.%02d", 5*i, 5*(i+1)));
        h_mass->GetXaxis()->SetTitle("Mass (GeV/c^{2})"); h_mass->GetYaxis()->SetTitle("Acceptance");
        h_mass->GetXaxis()->CenterTitle(true); h_mass->GetYaxis()->CenterTitle(true);
        h_mass->Draw("E1"); h_shivangi->Draw("E1 SAME");
        gPad->Update();

        std::vector<double> leg_coords = findEmptyRectangle(h_mass, h_shivangi);
        //TLegend *leg = new TLegend(leg_coords[0], leg_coords[1], leg_coords[2], leg_coords[3]);
        TLegend *leg = new TLegend(0.4, leg_coords[1], 0.8, leg_coords[3]);
        leg->SetBorderSize(0);
        leg->AddEntry(h_mass, "Latest Acceptance", "l");
        leg->AddEntry(h_shivangi, "Shivangi's Acceptance", "l");
        leg->Draw();
        c->SaveAs(Form("LD2_acceptance_xF_bin_%d.pdf", i));
        delete c;
    }
    f_mass->Close(); f_d_shivangi->Close();
}

void plot_combined_acceptance() {
    TFile *f_mass = TFile::Open("acceptance_mass_xF.root");
    TFile *f_h_shivangi = TFile::Open("acceptance_h_shivangi.root");
    TFile *f_d_shivangi = TFile::Open("acceptance_d_shivangi.root");
    if (!f_mass || !f_h_shivangi || !f_d_shivangi) return;

    for (int i = 0; i < 16; ++i) {
        TCanvas *c = new TCanvas(Form("c_comb_%d", i), Form("Combined Acceptance xF bin %d", i), 800, 600);
        c->SetLeftMargin(0.15);
        c->SetTickx(1); c->SetTicky(1);
        TH1F *h_mass = (TH1F*)f_mass->Get(Form("h_ratio_combine_xF_bin%d", i));
        TH1F *h_lh2 = (TH1F*)f_h_shivangi->Get(Form("LH2_acc_M_xF_%d", i));
        TH1F *h_ld2 = (TH1F*)f_d_shivangi->Get(Form("LD2_acc_M_xF_%d", i));
        if (!h_mass || !h_lh2 || !h_ld2) continue;

        TH1F *h_comb = (TH1F*)h_lh2->Clone(Form("h_comb_shivangi_%d", i));
        h_comb->Add(h_ld2); h_comb->Scale(0.5);
        h_mass->SetStats(0); h_comb->SetStats(0);
        h_mass->SetLineColor(kRed); h_comb->SetLineColor(kBlue);
        double maxY = 0;
        for (int bin = 1; bin <= h_mass->GetNbinsX(); ++bin) {
             maxY = std::max({maxY, h_mass->GetBinContent(bin) + h_mass->GetBinError(bin), h_comb->GetBinContent(bin) + h_comb->GetBinError(bin)});
        }
        h_mass->GetYaxis()->SetRangeUser(0, maxY * 1.1);
        h_mass->SetTitle(Form("Combined Acceptance for 0.%02d #leq xF < 0.%02d", 5*i, 5*(i+1)));
        h_mass->GetXaxis()->SetTitle("Mass (GeV/c^{2})"); h_mass->GetYaxis()->SetTitle("Acceptance");
        h_mass->GetXaxis()->CenterTitle(true); h_mass->GetYaxis()->CenterTitle(true);
        h_mass->Draw("E1"); h_comb->Draw("E1 SAME");
        gPad->Update();

        std::vector<double> leg_coords = findEmptyRectangle(h_mass, h_comb);
        //TLegend *leg = new TLegend(leg_coords[0], leg_coords[1], leg_coords[2], leg_coords[3]);
        TLegend *leg = new TLegend(0.4, leg_coords[1], 0.8, leg_coords[3]);
        leg->SetBorderSize(0);
        leg->AddEntry(h_mass, "Latest Combined Acceptance", "l");
        leg->AddEntry(h_comb, "Shivangi's Combined Acceptance", "l");
        leg->Draw();
        c->SaveAs(Form("Combined_acceptance_xF_bin_%d.pdf", i));
        delete c;
    }
    f_mass->Close(); f_h_shivangi->Close(); f_d_shivangi->Close();
}

void plot_acceptance_ratio() {
    TFile *f_mass = TFile::Open("acceptance_mass_xF.root");
    TFile *f_h_shivangi = TFile::Open("acceptance_h_shivangi.root");
    TFile *f_d_shivangi = TFile::Open("acceptance_d_shivangi.root");
    if (!f_mass || !f_h_shivangi || !f_d_shivangi) return;

    for (int i = 0; i < 16; ++i) {
        TCanvas *c = new TCanvas(Form("c_ratio_%d", i), Form("Acceptance Ratio xF bin %d", i), 800, 600);
        c->SetLeftMargin(0.15);
        c->SetTickx(1); c->SetTicky(1);
        TH1F *h_mass_r = (TH1F*)f_mass->Get(Form("h_ratio_acceptance_xF_bin%d", i));
        TH1F *h_lh2 = (TH1F*)f_h_shivangi->Get(Form("LH2_acc_M_xF_%d", i));
        TH1F *h_ld2 = (TH1F*)f_d_shivangi->Get(Form("LD2_acc_M_xF_%d", i));
        if (!h_mass_r || !h_lh2 || !h_ld2) continue;

        TH1F *h_shivangi_r = (TH1F*)h_lh2->Clone(Form("h_shivangi_ratio_%d", i));
        h_shivangi_r->Divide(h_ld2);
        h_mass_r->SetStats(0); h_shivangi_r->SetStats(0);
        h_mass_r->SetLineColor(kRed); h_shivangi_r->SetLineColor(kBlue);
        
        double maxY = -1e9, minY = 1e9;
        for (int bin = 1; bin <= h_mass_r->GetNbinsX(); ++bin) {
            maxY = std::max({maxY, h_mass_r->GetBinContent(bin) + h_mass_r->GetBinError(bin), h_shivangi_r->GetBinContent(bin) + h_shivangi_r->GetBinError(bin)});
            minY = std::min({minY, h_mass_r->GetBinContent(bin) - h_mass_r->GetBinError(bin), h_shivangi_r->GetBinContent(bin) - h_shivangi_r->GetBinError(bin)});
        }
        double range = maxY - minY; if (range <= 0) range = 1;
        h_mass_r->GetYaxis()->SetRangeUser(minY - 0.1 * range, maxY + 0.1 * range);

        h_mass_r->SetTitle(Form("Acceptance Ratio (LH2/LD2) for 0.%02d #leq xF < 0.%02d", 5*i, 5*(i+1)));
        h_mass_r->GetXaxis()->SetTitle("Mass (GeV/c^{2})"); h_mass_r->GetYaxis()->SetTitle("Ratio");
        h_mass_r->GetXaxis()->CenterTitle(true); h_mass_r->GetYaxis()->CenterTitle(true);
        h_mass_r->Draw("E1"); h_shivangi_r->Draw("E1 SAME");
        gPad->Update();

        std::vector<double> leg_coords = findEmptyRectangle(h_mass_r, h_shivangi_r);
        //TLegend *leg = new TLegend(leg_coords[0], leg_coords[1], leg_coords[2], leg_coords[3]);
        TLegend *leg = new TLegend(0.4, leg_coords[1], 0.8, leg_coords[3]);
        leg->SetBorderSize(0);
        leg->AddEntry(h_mass_r, "Latest Ratio (LH2/LD2)", "l");
        leg->AddEntry(h_shivangi_r, "Shivangi's Ratio (LH2/LD2)", "l");
        leg->Draw();
        c->SaveAs(Form("Acceptance_ratio_xF_bin_%d.pdf", i));
        delete c;
    }
    f_mass->Close(); f_h_shivangi->Close(); f_d_shivangi->Close();
}

void run_all_plots() {
    plot_LH2_acceptance();
    plot_LD2_acceptance();
    plot_combined_acceptance();
    plot_acceptance_ratio();
}
