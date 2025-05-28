import datetime
import subprocess
import os
import shutil # For checking if pdflatex is available

# --- Configuration ---
AUTHOR_NAME = "Weekly Updates"
UNIVERSITY_NAME = "New Mexico State University"
MAIN_TITLE_TEXT = "MEASUREMENT OF DY ABSOLUTE CROSS-SECTION MEASUREMENT FOR P+P AND P+D COLLISIONS WITH 120GeV PROTON BEAM AT FERMILAB"

# Define sections for the Table of Contents and slides
SECTIONS = [
    {
        "title": "Overview",
        "bullets": [
            "Event Selection: Chuck cuts.",
            "Files Used:",
            {"raw_latex_block": r"""% Note: This block is intended to be a sub-list or detailed text for "Files Used:"
                \begin{itemize}
                    \itemsep0.25em % Optional: add some space between sub-items
                    \item Data:\\ \texttt{roadset57\_70\_R008\_2111v42\_tmp\_noPhys.root}
                    \item $J/\psi$ MC Events:\\ \texttt{mc\_jpsi\_LH2\_M027\_S001\_messy\_occ\_pTxFweight\_v2.root}
                    \item $\psi'$ MC Events: \texttt{mc\_psiprime\_LH2\_M027\_S001\_messy\_occ\_pTxFweight\_v2.root}
                    \item Mixed Events:\\ \texttt{mixFPGA4\_67\_R008\_preBin2.root}
                    \item NNPDF4 File:\\ \texttt{NNPDF40\_xFnew\_p.root}
                    \item CT18 File:\\ \texttt{CT18\_xFnew\_p.root}
                \end{itemize}"""
            }
        ]
    },
    {
        "title": "Analysis Steps",
        "bullets": [
            "Data selection and quality checks.",
            "Background subtraction techniques.",
            "Acceptance and efficiency corrections.",
            "Systematic uncertainty estimation."
        ]
    },
    # New section for the PDF plots
    {
        "title": "Differential Cross-Section Plots", # This title is for the ToC section
        "is_plot_section": True # Custom flag to identify this section
    }
]

LOGO_FILENAME = "NMSU_logo.png" # Place this in the same directory as the script
PDF_PLOT_DIR = "../../src/Archived/DiffCross/plots/"

# xF Bin Definitions
XF_BINS = [
    (0.0, 0.05), (0.05, 0.1), (0.1, 0.15), (0.15, 0.2),
    (0.2, 0.25), (0.25, 0.3), (0.3, 0.35), (0.35, 0.4),
    (0.4, 0.45), (0.45, 0.5), (0.5, 0.55), (0.55, 0.6),
    (0.6, 0.65), (0.65, 0.7), (0.7, 0.75), (0.75, 0.8),
    (0.8, 0.85)
]


def get_formatted_date():
    """Returns the current date formatted as YYYY-MM-DD."""
    return datetime.datetime.now().strftime("%Y-%m-%d")

def tex_escape(text):
    """Escapes special LaTeX characters in a string meant for normal text output."""
    if not isinstance(text, str):
        try:
            text = str(text)
        except:
            return text # Or raise an error if non-convertible types are problematic
    return text.replace('&', '\\&').replace('%', '\\%').replace('$', '\\$') \
        .replace('#', '\\#').replace('_', '\\_').replace('{', '\\{') \
        .replace('}', '\\}').replace('~', '\\textasciitilde{}') \
        .replace('^', '\\textasciicircum{}').replace('\\', '\\textbackslash{}')

def generate_latex_code(author, institute, title, date_str, sections_data, logo_filename, plot_dir, xf_bins_data):
    """Generates the LaTeX Beamer code for the presentation."""

    author_esc = tex_escape(author)
    institute_esc = tex_escape(institute)
    title_esc = tex_escape(title) # Main title should be escaped

    footer_box_ht = "4.5ex"
    footer_box_dp = "2.5ex"
    footer_edge_padding = "1.2cm"
    footer_inter_padding = "0.3cm"

    latex_code = f"""
\\documentclass[10pt]{{beamer}} % Smaller base font size
\\usepackage[utf8]{{inputenc}}
\\usepackage{{lmodern}} % For better font rendering
\\usepackage{{textcomp}} % For text symbols
\\usepackage{{graphicx}} % For including images/PDFs
\\usepackage{{hyperref}} % For clickable links
\\usepackage{{amsmath}} % For math environments if needed
\\usepackage{{amsfonts}} % For math fonts if needed
\\usepackage{{amssymb}} % For math symbols if needed
\\usepackage{{xcolor}} % Required for \definecolor
\\usepackage{{ragged2e}} % For \justify command if needed

% --- NMSU Color Definitions ---
\\definecolor{{NMSUCrimson}}{{RGB}}{{140,11,66}}
\\definecolor{{NMSUWhite}}{{RGB}}{{255,255,255}}
\\definecolor{{NMSUGray}}{{RGB}}{{167,168,170}}
\\definecolor{{NMSUDarkGray}}{{RGB}}{{80,80,80}}

% --- Beamer Theme & Color Customization for NMSU ---
\\usetheme{{Madrid}}
\\setbeamertemplate{{navigation symbols}}{{}} % Remove navigation buttons/symbols

% Apply NMSU Colors
\\setbeamercolor{{normal text}}{{fg=NMSUDarkGray,bg=NMSUWhite}}
\\setbeamercolor{{structure}}{{fg=NMSUCrimson}}
\\setbeamercolor{{title}}{{fg=NMSUCrimson,bg=NMSUWhite}}
\\setbeamercolor{{author}}{{fg=NMSUDarkGray,bg=NMSUWhite}}
\\setbeamercolor{{institute}}{{fg=NMSUDarkGray,bg=NMSUWhite}}
\\setbeamercolor{{date}}{{fg=NMSUDarkGray,bg=NMSUWhite}}
\\setbeamercolor{{frametitle}}{{fg=NMSUWhite,bg=NMSUCrimson}}
\\setbeamercolor{{block title}}{{fg=NMSUWhite,bg=NMSUCrimson}}
\\setbeamercolor{{block body}}{{bg=NMSUGray!20}}
\\setbeamercolor{{alerted text}}{{fg=NMSUCrimson}}
\\setbeamercolor{{section in toc}}{{fg=NMSUDarkGray}}
\\setbeamercolor{{subsection in toc}}{{fg=NMSUDarkGray}}

% Font sizes
\\setbeamerfont{{title}}{{size=\\large}}
\\setbeamerfont{{author}}{{size=\\normalsize}}
\\setbeamerfont{{institute}}{{size=\\normalsize}}
\\setbeamerfont{{frametitle}}{{size=\\normalsize}}

% --- Presentation Metadata ---
\\title{{{title_esc}}}
\\author{{{author_esc}}}
\\institute{{{institute_esc}}}
\\date{{{date_str}}}

% --- Custom Title Page Template ---
\\setbeamertemplate{{title page}}{{
  \\vbox{{}}
  \\vfill
  \\begingroup
    \\centering
    \\begin{{beamercolorbox}}[sep=8pt,center]{{title}}
      \\usebeamerfont{{title}}\\inserttitle\\par%
      \\ifx\\insertsubtitle\\@empty\\else%
        \\vskip0.25em%
        {{\\usebeamerfont{{subtitle}}\\usebeamercolor[fg]{{subtitle}}\\insertsubtitle\\par}}%
      \\fi%
    \\end{{beamercolorbox}}%
    \\vskip1em\\par
    \\begin{{beamercolorbox}}[sep=8pt,center]{{author}}
      \\usebeamerfont{{author}}\\insertauthor
    \\end{{beamercolorbox}}
    \\vskip1.0em
    \\begin{{center}}
      \\includegraphics[height=2cm]{{{logo_filename}}}
    \\end{{center}}
  \\endgroup
  \\vfill
}}

% --- Custom Footline ---
\\setbeamercolor{{footline}}{{fg=NMSUWhite,bg=NMSUCrimson}}
\\setbeamercolor{{author in head/foot}}{{parent=footline}}
\\setbeamercolor{{institute in head/foot}}{{parent=footline}}
\\setbeamercolor{{date in head/foot}}{{parent=footline}}

\\setbeamertemplate{{footline}}{{
  \\leavevmode%
  \\hbox{{%
  \\begin{{beamercolorbox}}[wd=.333333\\paperwidth,ht={footer_box_ht},dp={footer_box_dp},leftskip={footer_edge_padding},rightskip={footer_inter_padding}]{{author in head/foot}}%
    \\usebeamerfont{{author in head/foot}}{{{author_esc}}}
  \\end{{beamercolorbox}}%
  \\begin{{beamercolorbox}}[wd=.333333\\paperwidth,ht={footer_box_ht},dp={footer_box_dp},leftskip={footer_inter_padding},rightskip={footer_inter_padding},center]{{institute in head/foot}}%
    \\usebeamerfont{{institute in head/foot}}{{{institute_esc}}}
  \\end{{beamercolorbox}}%
  \\begin{{beamercolorbox}}[wd=.333333\\paperwidth,ht={footer_box_ht},dp={footer_box_dp},leftskip={footer_inter_padding},rightskip={footer_edge_padding}]{{date in head/foot}}%
    \\hfill\\usebeamerfont{{date in head/foot}}{{{date_str}}}\\ifnum\\value{{framenumber}}>1 ~ \\insertframenumber\\,/\\,\\inserttotalframenumber\\fi
  \\end{{beamercolorbox}}%
  }}%
  \\vskip0pt%
}}
\\setbeamerfont{{author in head/foot}}{{size=\\scriptsize}}
\\setbeamerfont{{institute in head/foot}}{{size=\\scriptsize}}
\\setbeamerfont{{date in head/foot}}{{size=\\scriptsize}}

\\begin{{document}}

% --- Title Frame ---
\\begin{{frame}}
  \\titlepage
\\end{{frame}}

% --- Table of Contents Frame ---
\\begin{{frame}}{{Table of Contents}}
  \\tableofcontents
\\end{{frame}}
"""
    # --- Generate Section Frames ---
    for section_info in sections_data:
        # Section titles for ToC should be escaped if they might contain special chars
        section_title_for_toc = tex_escape(section_info["title"])
        latex_code += f"\n% --- {section_title_for_toc} Section ---\n"
        latex_code += f"\\section{{{section_title_for_toc}}}\n\n"

        if section_info.get("is_plot_section"):
            for i in range(16): # For LH2_0_cs.pdf to LH2_15_cs.pdf
                plot_filename = f"LH2_{i}_cs.pdf"
                relative_plot_path = os.path.join(plot_dir, plot_filename).replace("\\", "/")
                
                # Construct the frame title for plot slides (DO NOT tex_escape this specific string)
                frame_title_for_plot = f"Plot: LH2 {i} cs" # Default
                if 0 <= i < len(xf_bins_data):
                    lower_xf, upper_xf = xf_bins_data[i]
                    # This string is already valid LaTeX for the title
                    frame_title_for_plot = f"Double Differential Cross-Section: ${lower_xf} \\le x_F < {upper_xf}$"
                else:
                    print(f"Warning: xF bin number {i} is out of range for XF_BINS data. Using default title for {plot_filename}")
                
                latex_code += f"\\begin{{frame}}{{{frame_title_for_plot}}}\n" # Use directly
                latex_code += f"  \\centering\n"
                latex_code += f"  \\includegraphics[width=0.95\\textwidth, height=0.8\\textheight, keepaspectratio]{{{relative_plot_path}}}\n"
                latex_code += f"\\end{{frame}}\n\n"
        else:
            # For regular section titles on frames, they should also be escaped
            frame_title_on_slide = tex_escape(section_info["title"])
            latex_code += f"\\begin{{frame}}{{{frame_title_on_slide}}}\n"
            # latex_code += f"  \\frametitle{{{frame_title_on_slide}}}\n" # frametitle is already part of \begin{frame}
            
            bullets = section_info.get("bullets", [])
            if bullets:
                latex_code += f"  \\begin{{itemize}}\n"
                latex_code += f"    \\itemsep0.5em\n"
                for bullet_item in bullets:
                    if isinstance(bullet_item, str):
                        latex_code += f"    \\item {tex_escape(bullet_item)}\n"
                    elif isinstance(bullet_item, dict) and "raw_latex_block" in bullet_item:
                        latex_code += f"    {bullet_item['raw_latex_block']}\n"
                    else:
                        print(f"Warning: Unknown bullet type in section '{section_title_for_toc}': {type(bullet_item)}")
                latex_code += f"  \\end{{itemize}}\n"
            else:
                latex_code += f"  % No bullets for this section.\n"
            latex_code += f"\\end{{frame}}\n"

    latex_code += """
\\end{document}
"""
    return latex_code

def compile_latex_to_pdf(tex_filename):
    """Compiles the .tex file to .pdf using pdflatex."""

    if not shutil.which("pdflatex"):
        print("Error: 'pdflatex' command not found.")
        print("Please ensure a LaTeX distribution (TeX Live, MiKTeX, MacTeX) is installed and in your system PATH.")
        return None

    tex_dir = os.path.dirname(os.path.abspath(tex_filename))
    tex_file_basename = os.path.basename(tex_filename)

    original_dir = os.getcwd()
    os.chdir(tex_dir)

    compilation_success = False
    for i in range(2):
        try:
            process = subprocess.run(
                ['pdflatex', '-interaction=nonstopmode', '-file-line-error', tex_file_basename],
                capture_output=True, text=True, encoding='utf-8', check=False
            )
            if process.returncode != 0:
                print(f"Error during pdflatex run {i+1} (return code: {process.returncode}):")
                stdout_summary = process.stdout.strip()
                stderr_summary = process.stderr.strip()
                if stdout_summary:
                     print("--- STDOUT (summary) ---")
                     print((stdout_summary[:1000] + '...') if len(stdout_summary) > 1000 else stdout_summary)
                if stderr_summary:
                    print("--- STDERR (summary) ---")
                    print((stderr_summary[:1000] + '...') if len(stderr_summary) > 1000 else stderr_summary)
                compilation_success = False
                break
            else:
                print(f"pdflatex run {i+1} successful.")
                compilation_success = True
        except FileNotFoundError:
            print("Error: 'pdflatex' command not found. Please ensure LaTeX is installed and in your PATH.")
            os.chdir(original_dir)
            return None
        except Exception as e:
            print(f"An unexpected error occurred during pdflatex execution: {e}")
            os.chdir(original_dir)
            return None

    final_pdf_path = os.path.abspath(f"{os.path.splitext(tex_file_basename)[0]}.pdf")

    if compilation_success:
        print(f"PDF generated: {final_pdf_path}")
        extensions_to_clean = ['.aux', '.log', '.nav', '.snm', '.toc', '.out', '.lof', '.lol', '.lot', '.fls', '.fdb_latexmk']
        for ext in extensions_to_clean:
            aux_file = f"{os.path.splitext(tex_file_basename)[0]}{ext}"
            if os.path.exists(aux_file):
                try:
                    os.remove(aux_file)
                except OSError as e_clean:
                    print(f"Warning: Could not remove auxiliary file {aux_file}: {e_clean}")
        os.chdir(original_dir)
        return final_pdf_path
    else:
        print(f"PDF generation failed. Check the output from pdflatex above.")
        print(f"The .tex file '{os.path.abspath(tex_file_basename)}' and log files are available for inspection.")
        os.chdir(original_dir)
        return None

def main():
    current_date = get_formatted_date()
    output_filename_base = f"DY_CrossSection_NMSU_WeeklyUpdate_{current_date}"
    tex_filename = f"{output_filename_base}.tex"

    print(f"Generating LaTeX code for {tex_filename} with PDF plots...")
    latex_content = generate_latex_code(
        author=AUTHOR_NAME,
        institute=UNIVERSITY_NAME,
        title=MAIN_TITLE_TEXT,
        date_str=current_date,
        sections_data=SECTIONS,
        logo_filename=LOGO_FILENAME,
        plot_dir=PDF_PLOT_DIR,
        xf_bins_data=XF_BINS # Pass the XF_BINS data
    )

    try:
        with open(tex_filename, "w", encoding="utf-8") as f:
            f.write(latex_content)
        abs_tex_filename = os.path.abspath(tex_filename)
        print(f"LaTeX file '{abs_tex_filename}' created successfully.")
        print(f"IMPORTANT:")
        print(f"  Logo ('{LOGO_FILENAME}') should be in: {os.path.dirname(abs_tex_filename)}")
        print(f"  Plots should be in a directory accessible via the relative path '{PDF_PLOT_DIR}' from: {os.path.dirname(abs_tex_filename)}")
        example_plot_resolved_path = os.path.normpath(os.path.join(os.path.dirname(abs_tex_filename), PDF_PLOT_DIR, "LH2_0_cs.pdf"))
        print(f"  e.g., LaTeX will look for the first plot at: {example_plot_resolved_path}")

    except IOError as e:
        print(f"Error writing .tex file: {e}")
        return

    print("\nAttempting to compile LaTeX to PDF...")
    pdf_file_path = compile_latex_to_pdf(tex_filename)

    if pdf_file_path:
        print(f"\nSuccessfully generated PDF: '{pdf_file_path}'")
    else:
        abs_tex_filename_on_fail = os.path.abspath(tex_filename)
        print(f"\nFailed to generate PDF from '{abs_tex_filename_on_fail}'.")
        print("Please check the console output for errors from 'pdflatex'.")
        print("Ensure you have a working LaTeX distribution installed, the logo file is present, and the PDF plot paths are correct.")

if __name__ == '__main__':
    main()