# **Contributing to the Particle-In-Cell Scalable Application Resource (PICSAR)**


Please, read carefully the guidelines below when editing Fortran Files and 
contributing to the PXR repository. Respecting these guidelines ensures better readability of
the code by reviewers and future developpers. Below, Fortran blocks refer to portions of the code 
delimited by opening and closing statements. This includes control statements (IF/ELSE/END IF, 
DO/END DO, DO WHILE/END DO WHILE etc.) and Fortran constructs (MODULE/END MODULE, TYPE/END TYPE, 
SUBROUTINE/END SUBROUTINE, FUNCTION/END FUNCTION etc.)

## **Guidelines**

- Indentation width  of Fortran blocks has to be 2 white spaces (you can set this in your code editor),
  The parser in `utils/fortran_parsers/indent_file.py` can be used to automatically fixed this.  

- Fortran intrinsics and key words have to be written in UPPER CASE (e.g INTEGER, REAL, 
  SUBROUTINE, END, FUNCTION etc.). All other variables and custom names 
  must be written in lower case, 

- When creating  a fortran block of type **BLOCKTYPE* (e.g BLOCKTYPE can be MODULE, SUBROUTINE etc.) and name **blockname**, always end the block using 
"END BLOCKTYPE blockname" for readability,

- When using commas, always insert a space after the comma

- Lines cannot exceed 85 characters and ampersand continuation symbol should be added 
  at column 86 when continuing a line. The parser in `utils/fortran_parsers/justify_file.py` 
  can be used to automatically fixed this. 

- Keep a blank line between declaration blocks and code blocks for clarity 

