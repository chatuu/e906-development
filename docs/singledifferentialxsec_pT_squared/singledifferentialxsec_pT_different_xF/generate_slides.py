import os

# Configuration
BASE_DIR = "/root/github/e906-development/src/xsec_pT_squard/RS57-70"
OUTPUT_TEX = "slides.tex"

plot_files = [
    "Combined_XSec_Ratio_geom_pT.pdf",
    "Combined_XSec_Ratio_geom_pT2.pdf",
    "Combined_XSec_Ratio_true_pt_pT.pdf",
    "Combined_XSec_Ratio_true_pt_pT2.pdf",
    "CrossSection_LD2_geom_pT2_with_logo.pdf",
    "CrossSection_LD2_geom_pT_with_logo.pdf",
    "CrossSection_LD2_true_pt_pT2_with_logo.pdf",
    "CrossSection_LD2_true_pt_pT_with_logo.pdf",
    "CrossSection_LH2_geom_pT2_with_logo.pdf",
    "CrossSection_LH2_geom_pT_with_logo.pdf",
    "CrossSection_LH2_true_pt_pT2_with_logo.pdf",
    "CrossSection_LH2_true_pt_pT_with_logo.pdf",
    "CrossSection_Ratio_pd_2pp_geom_pT2_with_logo.pdf",
    "CrossSection_Ratio_pd_2pp_geom_pT_with_logo.pdf",
    "CrossSection_Ratio_pd_2pp_true_pt_pT2_with_logo.pdf",
    "CrossSection_Ratio_pd_2pp_true_pt_pT_with_logo.pdf"
]

def generate_tex():
    with open(OUTPUT_TEX, "w") as f:
        # Preamble
        f.write(r"""\documentclass{beamer}
\usetheme{Madrid}
\usepackage{graphicx}

\title{Addressing Comments Received}
\subtitle{$\frac{d\sigma}{d p^{2}_{T}}$ plots}
\author{Chatura Kuruppu}
\date{\today}

\begin{document}

\begin{frame}
    \titlepage
\end{frame}
""")

        # Generate Frames
        for filename in plot_files:
            # Clean title from filename
            clean_title = filename.replace(".pdf", "").replace("_with_logo", "").replace("_", " ")
            
            f.write(r"""
\begin{frame}{""" + clean_title + r"""}
    \centering
    \includegraphics[width=0.9\textwidth,height=0.7\textheight,keepaspectratio]{""" + os.path.join(BASE_DIR, filename) + r"""}
\end{frame}
""")

        f.write(r"\end{document}")

if __name__ == "__main__":
    generate_tex()
    print(f"Successfully generated {OUTPUT_TEX}")