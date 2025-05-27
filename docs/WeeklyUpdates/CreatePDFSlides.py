import datetime
import subprocess
import os
import shutil # For checking if pdflatex is available

# --- Configuration ---
AUTHOR_NAME = "Weekly Updates" 
UNIVERSITY_NAME = "New Mexico State University"
MAIN_TITLE_TEXT = "MEASUREMENT OF DY ABSOLUTE CROSS-SECTION MEASUREMENT FOR P+P AND P+D COLLISIONS WITH 120 GEV PROTON BEAM AT FERMILAB"
OVERVIEW_TITLE_TEXT = "Overview"
OVERVIEW_BULLETS = [
    "This is the first bullet point for the overview.",
    "This is the second bullet point, providing more details.",
    "And finally, the third bullet point summarizing key aspects."
]
LOGO_FILENAME = "NMSU_logo.png" # <<< Make sure this logo file exists in the script's directory

def get_formatted_date():
    """Returns the current date formatted as<y_bin_46>-MM-DD."""
    return datetime.datetime.now().strftime("%Y-%m-%d")

def tex_escape(text):
    """Escapes special LaTeX characters in a string."""
    return text.replace('&', '\\&').replace('%', '\\%').replace('$', '\\$') \
               .replace('#', '\\#').replace('_', '\\_').replace('{', '\\{') \
               .replace('}', '\\}').replace('~', '\\textasciitilde{}') \
               .replace('^', '\\textasciicircum{}').replace('\\', '\\textbackslash{}')

def generate_latex_code(author, institute, title, date_str, overview_title, bullets, logo_filename):
    """Generates the LaTeX Beamer code for the presentation."""
    
    author_esc = tex_escape(author)
    institute_esc = tex_escape(institute)
    title_esc = tex_escape(title)
    overview_title_esc = tex_escape(overview_title)
    bullets_esc = [tex_escape(b) for b in bullets]

    footer_box_ht = "4.5ex"
    footer_box_dp = "2.5ex"
    footer_edge_padding = "1.2cm"
    footer_inter_padding = "0.3cm"

    latex_code = f"""
\\documentclass{{beamer}}
\\usepackage[utf8]{{inputenc}}
\\usepackage{{lmodern}}
\\usepackage{{textcomp}}
\\usepackage{{graphicx}}

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

% --- Presentation Metadata ---
\\title{{{title_esc}}}
\\author{{{author_esc}}} 
\\institute{{{institute_esc}}}
\\date{{{date_str}}}
% DO NOT define \\logo globally if it should only appear on the title page via custom template

% --- Custom Title Page Template ---
% This redefines how \\titlepage looks.
% It includes Title, Author, and directly embeds the logo image,
% but omits Institute and Date from the main body.
% MODIFIED: Removed shadow=true, rounded=true from title and author beamercolorboxes
\\setbeamertemplate{{title page}}{{
  \\vbox{{}}
  \\vfill
  \\begingroup
    \\centering
    % Title
    \\begin{{beamercolorbox}}[sep=8pt,center]{{title}} % Removed shadow & rounded options
      \\usebeamerfont{{title}}\\inserttitle\\par%
      \\ifx\\insertsubtitle\\@empty\\else%
        \\vskip0.25em%
        {{\\usebeamerfont{{subtitle}}\\usebeamercolor[fg]{{subtitle}}\\insertsubtitle\\par}}%
      \\fi%
    \\end{{beamercolorbox}}%
    \\vskip1em\\par
    % "Author" field (now "Weekly Updates")
    \\begin{{beamercolorbox}}[sep=8pt,center]{{author}} % Removed shadow & rounded options
      \\usebeamerfont{{author}}\\insertauthor
    \\end{{beamercolorbox}}
    % Direct Logo Placement
    \\vskip1.0em 
    \\begin{{center}}
      \\includegraphics[height=2.5cm]{{{logo_filename}}}
    \\end{{center}}
  \\endgroup
  \\vfill
}}

% --- Custom Footline: NMSU Colors, Thickness, Padding, Conditional Slide Numbers ---
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

% --- Title Frame (Frame Number 1) ---
\\begin{{frame}} 
  \\titlepage 
\\end{{frame}}

% --- Overview Frame (Frame Number 2 onwards) ---
\\begin{{frame}}
  \\frametitle{{{overview_title_esc}}}
  \\begin{{itemize}}
"""
    for item in bullets_esc:
        latex_code += f"    \\item {item}\n"
    latex_code += """
  \\end{itemize}
\\end{frame}

\\end{document}
"""
    return latex_code

def compile_latex_to_pdf(tex_filename):
    """Compiles the .tex file to .pdf using pdflatex."""
    
    if not shutil.which("pdflatex"):
        print("Error: 'pdflatex' command not found.")
        print("Please ensure a LaTeX distribution (TeX Live, MiKTeX, MacTeX) is installed and in your system PATH.")
        return None

    base_name = os.path.splitext(tex_filename)[0]
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
                stdout_summary = (process.stdout[:1000] + '...') if len(process.stdout) > 1000 else process.stdout
                stderr_summary = (process.stderr[:1000] + '...') if len(process.stderr) > 1000 else process.stderr
                print("--- STDOUT (summary) ---")
                print(stdout_summary)
                if stderr_summary: 
                    print("--- STDERR (summary) ---")
                    print(stderr_summary)
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
    
    pdf_filename_in_tex_dir = f"{os.path.splitext(tex_file_basename)[0]}.pdf"
    final_pdf_path = os.path.join(tex_dir, pdf_filename_in_tex_dir)

    if compilation_success:
        print(f"PDF generated: {final_pdf_path}")
        extensions_to_clean = ['.aux', '.log', '.nav', '.snm', '.toc', '.out', '.lof', '.lol', '.lot', '.fls', '.fdb_latexmk']
        for ext in extensions_to_clean:
            aux_file = f"{os.path.splitext(tex_file_basename)[0]}{ext}" 
            if os.path.exists(aux_file):
                try:
                    os.remove(aux_file)
                except OSError as e:
                    print(f"Warning: Could not remove auxiliary file {aux_file}: {e}")
        os.chdir(original_dir)
        return final_pdf_path 
    else:
        print(f"PDF generation failed. Check the output from pdflatex above.")
        print(f"The .tex file '{tex_file_basename}' and log files are available in '{tex_dir}' for inspection.")
        os.chdir(original_dir)
        return None

def main():
    current_date = get_formatted_date()
    
    output_filename_base = f"DY_CrossSection_NMSU_Latex_FlatTitle_{current_date}" # Updated filename
    tex_filename = f"{output_filename_base}.tex" 
    
    print(f"Generating LaTeX code for {tex_filename} with flat title page elements...")
    latex_content = generate_latex_code(
        author=AUTHOR_NAME,
        institute=UNIVERSITY_NAME,
        title=MAIN_TITLE_TEXT,
        date_str=current_date,
        overview_title=OVERVIEW_TITLE_TEXT,
        bullets=OVERVIEW_BULLETS,
        logo_filename=LOGO_FILENAME
    )
    
    try:
        with open(tex_filename, "w", encoding="utf-8") as f:
            f.write(latex_content)
        abs_tex_filename = os.path.abspath(tex_filename)
        print(f"LaTeX file '{abs_tex_filename}' created successfully.")
        print(f"IMPORTANT: Ensure your logo file ('{LOGO_FILENAME}') is in the same directory.")
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
        print("Ensure you have a working LaTeX distribution installed and the logo file is present.")

if __name__ == '__main__':
    main()