import fitz  # PyMuPDF
import sys
import os

def add_logo_to_pdf(pdf_input_path, logo_path, pdf_output_path):
    """
    Overlays a PNG logo onto the first page of a PDF file.
    """
    # <<< ISSUE 2: MODIFY THIS LINE to change the logo's position.
    # The format is fitz.Rect(left_x, top_y, right_x, bottom_y).
    # The page origin (0,0) is the TOP-LEFT corner.
    xpos = -360
    ypos = -169

    image_rectangle = fitz.Rect(300 + ypos, 200 + xpos, 360 + ypos, 315 + xpos)  # Adjusted position for the logo

    try:
        # Open the source PDF
        pdf_file = fitz.open(pdf_input_path)
        
        # Get the first page
        first_page = pdf_file[0]
        
        # <<< ISSUE 1: MODIFIED THIS LINE to rotate the logo by 90 degrees.
        first_page.insert_image(image_rectangle, filename=logo_path, rotate=90)
        
        # Save the result to the new output path
        pdf_file.save(pdf_output_path)
        pdf_file.close()
        
        print(f"✅ Successfully added logo to '{pdf_output_path}'")

    except Exception as e:
        print(f"❌ An error occurred: {e}")

# --- Main execution ---
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python add_logo.py <input_pdf_path> <logo_png_path>")
        sys.exit(1)

    input_pdf = sys.argv[1]
    logo_file = sys.argv[2]
    
    file_name, file_extension = os.path.splitext(input_pdf)
    output_pdf = f"{file_name}_with_logo{file_extension}"
    
    add_logo_to_pdf(input_pdf, logo_file, output_pdf)